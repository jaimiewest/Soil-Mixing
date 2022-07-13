### Load required packages
library(phyloseq)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(stringr)
library(doMC)
library(harrietr)
library(magrittr)
library(picante) # requires ape

#### Getting started #### 

# read in ps object
ps.TRT=readRDS("ps.CON18.remake")
# prune away taxa with 0 counts
ps.TRT.prune = prune_taxa(taxa_sums(ps.TRT) > 0, ps.TRT)

# normalize to relative abundances
ps.TRT.norm = transform_sample_counts(ps.TRT.prune, function(x) x / sum(x))
ps = ps.TRT.norm

# read in the phylogeny
phylo = ggtree::read.tree("tree_CON18_remake_Apr2021.nwk")
phy_tree(ps)=phylo

#### CALCULATE bMNTD for treatments #### 
# Goal here is to calculate bMNTD for every possible pair of samples in each mixing set
# Vortex controls are all compared to each other, within mixign frequency.

# Create empty data frame to hold all results
bMNTD.df = data.frame()
# Then, run loop for getting data for all TimesMixed except 1 and initial, without Vortex controls
for (i in c("2","4","8","16","32")){
  for (j in c(1,2,3,4)){
    # Get a phyloseq object with only our samples of interest
    ps.bMNTD = subset_samples(ps,sample_data(ps)$TimesMixed == i &
                                sample_data(ps)$Rep == j &
                                sample_data(ps)$VortexControl == "N")
    # Get OTU table
    otu_full=otu_table(ps.bMNTD)
    # Get all sample names
    sampleNames = sample_names(ps.bMNTD)
    # Create empty matrix
    resultsMatrix=matrix(rep(NA,length(sampleNames)^2), nrow=length(sampleNames), dimnames = list(c(sampleNames), c(sampleNames)))
    # Create all comparisons (each sample to every other sample)
    Comparisons = expand.grid(sampleNames,sampleNames)
    # Do multicore, if possible on your computer
    registerDoMC(cores=2)
    # Multicore populate results matrix
    resultsMatrix = foreach(i=1:dim(Comparisons)[1]) %dopar%{
      # Get the comparisons of interest
      com1 = Comparisons[i,1]
      com2 = Comparisons[i,2]
      # Drop otu table to only those communities
      otu = otu_full[c(com1,com2)]
      # Create phylogeny-otu table object from phylogenetic tree
      match.phylo.otu = picante::match.phylo.data(phylo, t(otu))
      # Make vector of names where no taxa are present in this set of sample comparisons
      tipstodrop = row.names(match.phylo.otu$data)[rowSums(match.phylo.otu$data)==0]
      # Dropping those tips
      match.phylo.otu$phy = drop.tip(match.phylo.otu$phy,tipstodrop)
      match.phylo.otu = picante::match.phylo.data(match.phylo.otu$phy,match.phylo.otu$data)
      
      # Calculate bMNTD; include conspecifics (i.e., if the same taxon is in both communities, dist is zero)
      beta.mntd.weighted = as.matrix(picante::comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T,exclude.conspecifics=F));
      beta.mntd.weighted
      resultsMatrix[com1,com2]=beta.mntd.weighted[2,1]
    }
    # Adjusting format of results matrix (probably inefficiently)
    resultsMatrixFinal= matrix(as.matrix(resultsMatrix,nrows=length(sampleNames)),ncol=length(sampleNames))
    colnames(resultsMatrixFinal)=sampleNames
    row.names(resultsMatrixFinal)=sampleNames
    resultsMatrixFinal = data.frame(resultsMatrixFinal)
    bMNTD = as.matrix(resultsMatrixFinal)
    row.names(bMNTD) = colnames(bMNTD)
    bMNTD = melt_dist(bMNTD, dist_name = "bMNTD")
    # Add columns with run info
    bMNTD$TimesMixed = i
    bMNTD$Rep = j
    bMNTD$VortexControl = "N"
    bMNTD.df = rbind(bMNTD.df,bMNTD)
  }
}
head(bMNTD.df)
dim(bMNTD.df) 
# Looks ok; this dataset should have 560 observations (28 possible combinations x 4 mixing sets x 5 treatments)

# Then, a very similar for loop to add on the vortex controls, which are compared to each other within each times mixed trt

for (i in c("2","4","8","16","32")){
  ps.bMNTD = subset_samples(ps,sample_data(ps)$TimesMixed == i &
                              sample_data(ps)$VortexControl == "Y")
  otu_full=otu_table(ps.bMNTD)
  sampleNames = sample_names(ps.bMNTD)
  resultsMatrix=matrix(rep(NA,length(sampleNames)^2), nrow=length(sampleNames), dimnames = list(c(sampleNames), c(sampleNames)))
  Comparisons = expand.grid(sampleNames,sampleNames)
  registerDoMC(cores=2)
  resultsMatrix = foreach(i=1:dim(Comparisons)[1]) %dopar%{
    com1 = Comparisons[i,1]
    com2 = Comparisons[i,2]
    otu = otu_full[c(com1,com2)]
    match.phylo.otu = picante::match.phylo.data(phylo, t(otu))
    tipstodrop = row.names(match.phylo.otu$data)[rowSums(match.phylo.otu$data)==0]
    match.phylo.otu$phy = drop.tip(match.phylo.otu$phy,tipstodrop)
    match.phylo.otu = picante::match.phylo.data(match.phylo.otu$phy,match.phylo.otu$data)
    
    # Calculate bMNTD
    beta.mntd.weighted = as.matrix(picante::comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T,exclude.conspecifics=F));
    beta.mntd.weighted
    resultsMatrix[com1,com2]=beta.mntd.weighted[2,1]
  }
  resultsMatrixFinal= matrix(as.matrix(resultsMatrix,nrows=length(sampleNames)),ncol=length(sampleNames))
  colnames(resultsMatrixFinal)=sampleNames
  row.names(resultsMatrixFinal)=sampleNames
  resultsMatrixFinal = data.frame(resultsMatrixFinal)
  bMNTD = as.matrix(resultsMatrixFinal)
  row.names(bMNTD) = colnames(bMNTD)
  bMNTD = melt_dist(bMNTD, dist_name = "bMNTD")
  bMNTD$TimesMixed = i
  bMNTD$Rep = NA
  bMNTD$VortexControl = "Y"
  bMNTD.df = rbind(bMNTD.df,bMNTD)
}
head(bMNTD.df)
dim(bMNTD.df) # now there should be 700 total (540 from trts, above, plus additional 28 Vortex comparisons x 5 trts)

# Adjust column names 
colnames(bMNTD.df)
colnames(bMNTD.df)[1:2] = c("com1","com2")

# Make sure all columns are in the right class
str(bMNTD.df)
bMNTD.df$bMNTD = as.numeric(bMNTD.df$bMNTD)
bMNTD.df$Rep = as.character(bMNTD.df$Rep)

#### CALCULATE bMNTD for null model ####
# Every possible comparison amongst 1X controls

bMNTD.null = data.frame()
ps.bMNTD = subset_samples(ps,sample_data(ps)$TimesMixed == "1" & 
                            sample_data(ps)$VortexControl == "N")
otu_full=otu_table(ps.bMNTD)
sampleNames = sample_names(ps.bMNTD)
resultsMatrix=matrix(rep(NA,length(sampleNames)^2), nrow=length(sampleNames), dimnames = list(c(sampleNames), c(sampleNames)))
Comparisons = expand.grid(sampleNames,sampleNames)
registerDoMC(cores=2)
resultsMatrix = foreach(i=1:dim(Comparisons)[1]) %dopar%{
  com1 = Comparisons[i,1]
  com2 = Comparisons[i,2]
  otu = otu_full[c(com1,com2)]
  match.phylo.otu = picante::match.phylo.data(phylo, t(otu))
  tipstodrop = row.names(match.phylo.otu$data)[rowSums(match.phylo.otu$data)==0]
  match.phylo.otu$phy = drop.tip(match.phylo.otu$phy,tipstodrop)
  match.phylo.otu = picante::match.phylo.data(match.phylo.otu$phy,match.phylo.otu$data)
  
  # Calculate bMNTD
  beta.mntd.weighted = as.matrix(picante::comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T,exclude.conspecifics=F));
  beta.mntd.weighted
  resultsMatrix[com1,com2]=beta.mntd.weighted[2,1]
}
resultsMatrixFinal= matrix(as.matrix(resultsMatrix,nrows=length(sampleNames)),ncol=length(sampleNames))
colnames(resultsMatrixFinal)=sampleNames
row.names(resultsMatrixFinal)=sampleNames
resultsMatrixFinal = data.frame(resultsMatrixFinal)
bMNTD = as.matrix(resultsMatrixFinal)
row.names(bMNTD) = colnames(bMNTD) 
bMNTD = melt_dist(bMNTD, dist_name = "bMNTD")
bMNTD$TimesMixed = "1"
bMNTD$Rep = NA
bMNTD$VortexControl = "N"
bMNTD.null = rbind(bMNTD.null,bMNTD)

head(bMNTD.null)
dim(bMNTD.null) #should be 496 total (all combinations of 32 samples in 1x)

# Adjust column names 
colnames(bMNTD.null)
colnames(bMNTD.null)[1:2] = c("com1","com2")

# Make sure all columns are in the right class
str(bMNTD.null)
bMNTD.null$bMNTD = as.numeric(bMNTD.null$bMNTD)
bMNTD.null$Rep = as.character(bMNTD.null$Rep)

# This should have produced: Two dataframes with all sample comparisons of interest
# For all comparisons within a mixing set, or within the same frequency of vortex mixing (for vortex controls)
# and all comparisons within 1X, for the null model.

#### CALCULATE Bray-Curtis dissimilarities ####

# Calculate bray-curtis dissimilarities between all samples (we don't need all comparisons, but it's quick)
BC = as.matrix(distance(ps,method="bray"))

# Make names match bMNTD df
row.names(BC) = paste("X",row.names(BC),sep="")
colnames(BC) = paste("X",colnames(BC),sep="")
BC = melt_dist(BC, dist_name = "Bray_Curtis")
colnames(BC)[1:2] = c("com1","com2")
str(BC)

#### MERGE THE BMNTD AND BC DISSIMILARILY DF'S ####
# Make a column listing the comparison, in each df
bMNTD.df$Comparison = paste(bMNTD.df$com1,bMNTD.df$com2)
bMNTD.null$Comparison = paste(bMNTD.null$com1,bMNTD.null$com2)
BC$Comparison = paste(BC$com1,BC$com2)

# Make sure every bMNTD comparison is in the BC dataframe
# (sum should equal first number in dim)
dim(bMNTD.df)
sum(bMNTD.df$Comparison %in% BC$Comparison)

dim(bMNTD.null)
sum(bMNTD.null$Comparison %in% BC$Comparison)

# and now, merge bMNTD and BC dataframes
bMNTD.BC = merge(bMNTD.df, BC, by.x="Comparison", by.y="Comparison", all.x=TRUE, all.y=FALSE)
bMNTD.BC.null = merge(bMNTD.null, BC, by.x="Comparison", by.y="Comparison", all.x=TRUE, all.y=FALSE)

# Check it out
tail(bMNTD.BC)
dim(bMNTD.BC)

tail(bMNTD.BC.null)
dim(bMNTD.BC.null)

### Clean up the dataframes
# Get rid of 'artifact' X's
bMNTD.BC$com1 = bMNTD.BC$com1.x
bMNTD.BC$com2 = bMNTD.BC$com2.x
bMNTD.BC$com1 <- gsub("^.{0,1}", "", bMNTD.BC$com1)
bMNTD.BC$com2 <- gsub("^.{0,1}", "", bMNTD.BC$com2)
bMNTD.BC$Comparison <- paste(bMNTD.BC$com1, bMNTD.BC$com2)

# Delete the .x and .y columns
# First, confirm which columns are com1.x, com2.x, com1.y, com2.y.
head(bMNTD.BC)
bMNTD.BC2 = bMNTD.BC[-c(2, 3, 8, 9)]
# confirm that .x and .y columns were dropped, and Bray_Curtis column is still there
head(bMNTD.BC2)
# All good, reset bMNTD.BC df
bMNTD.BC = bMNTD.BC2

# Do the same for the null dataset
bMNTD.BC.null$com1 = bMNTD.BC.null$com1.x
bMNTD.BC.null$com2 = bMNTD.BC.null$com2.x
bMNTD.BC.null$com1 <- gsub("^.{0,1}", "", bMNTD.BC.null$com1)
bMNTD.BC.null$com2 <- gsub("^.{0,1}", "", bMNTD.BC.null$com2)
bMNTD.BC.null$Comparison <- paste(bMNTD.BC.null$com1, bMNTD.BC.null$com2)

# Delete the .x and .y columns
# First, confirm which columns are com1.x, com2.x, com1.y, com2.y. (Do not count the row name column)
head(bMNTD.BC.null)
bMNTD.BC.null2 = bMNTD.BC.null[-c(2, 3, 8, 9)]
# confirm that .x and .y columns were dropped, and Bray_Curtis column is still there
head(bMNTD.BC.null2)
# All good, reset bMNTD.BC.null df
bMNTD.BC.null = bMNTD.BC.null2


# Make sure all columns are in the right class and order
str(bMNTD.BC)
##if not:
# bMNTD.BC$bMNTD = as.numeric(bMNTD.BC$bMNTD)
# bMNTD.BC$Bray_Curtis = as.numeric(bMNTD.BC$Bray_Curtis)
# bMNTD.BC$Rep = as.character(bMNTD.BC$Rep)
# bMNTD.BC$TimesMixed = as.character(bMNTD.BC$TimesMixed)
bMNTD.BC$TimesMixed = ordered(bMNTD.BC$TimesMixed,levels=c("1","2","4","8","16","32")) #for graphing later

str(bMNTD.BC.null)
# bMNTD.BC.null$bMNTD = as.numeric(bMNTD.BC.null$bMNTD)
# bMNTD.BC.null$Bray_Curtis = as.numeric(bMNTD.BC.null$Bray_Curtis)
# bMNTD.BC.null$TimesMixed = as.character(bMNTD.BC.null$TimesMixed)

## Option to save the dataframes
write.csv(bMNTD.BC,"Derived_data/bMNTD.Bray_Curtis.Within_Mixing_Set.csv")
write.csv(bMNTD.BC.null,"Derived_data/bMNTD.Bray_Curtis.Null.csv")


####  assign quantiles (for bMNTD and BC)

## Read in metadata
meta = read.csv("sample_metadata_connors_fall2018.2021.csv",
                header = TRUE)

# create variable for number of comparisons; this will be used to calculate quantiles
obs = nrow(bMNTD.BC.null)
                                      
# create a vector with even quantile steps
q = seq(0, ((obs-1)/obs), by=(1/obs))
# Assign quantiles (non-parametrically) to each comparison, based on bMNTD values
# first sort rows/observations by bMNTD values
inter.df = inter.df[order(inter.df$bMNTD),]
# then assign quantiles to ascending bMNTD
inter.df$bMNTD.quant = q
# Similarly, assign quantiles (non-paramentrically) to each comparison, based on Bray_Curtis values
inter.df = inter.df[order(inter.df$Bray_Curtis),]
inter.df$B_C.quant = q
  
#### Using the null set quantiles, assign quantiles for within-treatment (within mixing set) comparisons
  
# add bMNTD.quant and B_C.quant columns to bMNTD.BC dataframe
# which contains all of the within-mixing set comparisons for each treatment
bMNTD.BC$bMNTD.quant <- NA
bMNTD.BC$B_C.quant <- NA
  
# Create function that takes a number (from treatment df "bMNTD.BC"), and identifies the closest value 
# in the null set df "inter.df" and populates the trt df (bMNTD.BC) with corresponding quantile
closest <- function(reflist, trtvalue){
    reflist[which(abs(reflist-trtvalue)==min(abs(reflist-trtvalue)))]}
  
# find bMNTD quantiles
for(r in 1:(nrow(bMNTD.BC))) {
  bMNTD.BC$bMNTD.quant[r] = inter.df[inter.df$bMNTD==closest(inter.df$bMNTD[],bMNTD.BC$bMNTD[r]),]$bMNTD.quant
}
# find Bray_Curtis quantiles
for(r in 1:(nrow(bMNTD.BC))) {
  bMNTD.BC$B_C.quant[r] = inter.df[inter.df$Bray_Curtis==closest(inter.df$Bray_Curtis[],bMNTD.BC$Bray_Curtis[r]),]$B_C.quant
}
  
### Then, assign cutoff quantiles to community assembly processes. 
#create empty column for the process ID
bMNTD.BC$Process <- NA
# establish the lower and upper cutoffs, approx equivalent to 0.025th and 0.975th quantiles.
# We could set the cutoffs equal to the target quantile values, but depending on the number of comparisons,
# this won't be a round number; and, we have to account for the 1st quantile step = 0!
# so we will use "<=" for lower quantile, to capture the 0th quantile
# and ">" for upper quantile (not >=).
  
# Find the closest element of q to 0.05, i.e. how many elements get us closest to 0.05
howmany= which(abs(q-0.05)==min(abs(q-0.05)))
# Now make it an even number
howmany2 = if((howmany %% 2) == 0) {
  which(abs(q-0.05)==min(abs(q-0.05)))
} else {
  (which(abs(q-0.05)==min(abs(q-0.05))) - 1) # Could add 1 instead..
}
# Establish the lower quantile cutoff
lowerquant= q[howmany2/2]
# establish the upper quantile cutoff.
upperquant= q[nobs(q)-(howmany2/2)]
  
# assign process ID's, considering both bMNTD and BC-Dis cutoffs for each comparison
for(r in 1:(nrow(bMNTD.BC))) {
  if (bMNTD.BC$bMNTD.quant[r] <= lowerquant) {
    bMNTD.BC$Process[r] = "Homogeneous Selection"
  } else if (bMNTD.BC$B_C.quant[r] <= lowerquant && bMNTD.BC$bMNTD.quant[r] > lowerquant && bMNTD.BC$bMNTD.quant[r] <= upperquant) {
    bMNTD.BC$Process[r] = "Homogenizing Dispersal"
  } else if (bMNTD.BC$bMNTD.quant[r] > upperquant) {
    bMNTD.BC$Process[r] = "Variable Selection"
  } else if (bMNTD.BC$B_C.quant[r] > upperquant && bMNTD.BC$bMNTD.quant[r] > lowerquant && bMNTD.BC$bMNTD.quant[r] <= upperquant) {
    bMNTD.BC$Process[r] = "Dispersal Limitation"
  } else if (bMNTD.BC$B_C.quant[r] > lowerquant && bMNTD.BC$B_C.quant[r] <= upperquant && bMNTD.BC$bMNTD.quant[r] > lowerquant && bMNTD.BC$bMNTD.quant[r] <= upperquant) {
    bMNTD.BC$Process[r] = "Undominated"
  } else bMNTD.BC$Process[r] = NA
}
  
## Option to save data
write.csv(bMNTD.BC,"Derived_data/ProcessID_Assignments.csv")


#### TALLY/SUMMARIZE THE COMMUNITY ASSEMBLY PROCESSES ASSIGNED TO WITHIN-TREATMENT COMPARISONS ####

# Check that no ID assignments were missed; there should be no NAs for any of the Iter. columns
bMNTD.BC %>%
  select_if(function(x) any(is.na(x))) %>%
  summarise_each(funs(sum(is.na(.))))

# summarize the process ID assignments.
ProcessSummary = bMNTD.BC %>%
  group_by(TimesMixed,VortexControl,Process)%>%
  dplyr::summarize(Count = n())
head(ProcessSummary)
str(ProcessSummary)

#### STACKED BAR GRAPH--TREATMENT COMMUNITIES ####

# Rename dataframe, for ease of use
ProcID = data.frame(ProcessSummary)

### Subset for TREATMENT communities
ProcID.trt = subset(ProcID, VortexControl=="N")
ProcID.trt

# Order factors for graphing
ProcID.trt$TimesMixed = ordered(ProcID.trt$TimesMixed,levels=c("2", "4", "8", "16", "32"))
ProcID.trt$Process = ordered(ProcID.trt$Process, levels=c("Undominated",
                                            "Homogenizing Dispersal",
                                            "Homogeneous Selection",
                                            "Variable Selection",
                                            "Dispersal Limitation"))


### Stacked bar graph; COmmunity assembly processes within each treatment (by mixing set)
ID_Palette = c("#CCCCCC", "#336699", "#99CC99", "#FFCC00", "#FF9900", "#CC3300")
n.obs = 112
Trt.plot = ggplot(ProcID.trt, aes(x=factor(TimesMixed), y=Count/n.obs)) +
  geom_bar(position = "stack", stat="identity",aes(fill=Process)) +
  expand_limits(y=c(0.0,1.04)) +
  scale_fill_manual(values = ID_Palette)+
  scale_x_discrete(labels=c("2"="2\u00d7", "4"="4\u00d7", "8"="8\u00d7", "16"="16\u00d7", "32"="32\u00d7")) +
  labs(x = "Mixing frequency \nComparisons within pooled sets", y = "Proportion of comparisons", fill = "Community assembly process") +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+
  #theme(legend.position = "none") 
Trt.plot #600 x 500
            
            
            
            
            
            
            
#### Now to estimate the community assembly processes for each pooled sample + a vortex control combo, from the same mixing frequency
            
#### CALCULATE bMNTD for treatments #### 
# Goal here is to calculate bMNTD for every possible pair: a pooled sample + a vortex control of the same mixing frequency.

# Create empty data frame to hold all results
bMNTD.vs.vort = data.frame()
# Then, run loop for getting data for all TimesMixed except 1 and initial, without Vortex controls
for (i in c("2","4","8","16","32")){
    # Get a phyloseq object with only our samples of interest
    ps.bMNTD = subset_samples(ps,sample_data(ps)$TimesMixed == i)
    # Get OTU table
    otu_full=otu_table(ps.bMNTD)
    # Get all sample names
    sampleNames = sample_names(ps.bMNTD)
    # Create empty matrix
    resultsMatrix=matrix(rep(NA,length(sampleNames)^2), nrow=length(sampleNames), dimnames = list(c(sampleNames), c(sampleNames)))
    # Create all comparisons (each sample to every other sample)
    Comparisons = expand.grid(sampleNames,sampleNames)
    # Do multicore, if possible on your computer
    registerDoMC(cores=2)
    # Multicore populate results matrix
    resultsMatrix = foreach(i=1:dim(Comparisons)[1]) %dopar%{
      # Get the comparisons of interest
      com1 = Comparisons[i,1]
      com2 = Comparisons[i,2]
      # Drop otu table to only those communities
      otu = otu_full[c(com1,com2)]
      # Create phylogeny-otu table object from phylogenetic tree
      match.phylo.otu = picante::match.phylo.data(phylo, t(otu))
      # Make vector of names where no taxa are present in this set of sample comparisons
      tipstodrop = row.names(match.phylo.otu$data)[rowSums(match.phylo.otu$data)==0]
      # Dropping those tips
      match.phylo.otu$phy = drop.tip(match.phylo.otu$phy,tipstodrop)
      match.phylo.otu = picante::match.phylo.data(match.phylo.otu$phy,match.phylo.otu$data)
      
      # Calculate bMNTD... it is interesting to consider excluding conspecifics...
      # For now, we include conspecifics (i.e., if there is the same taxon in both communities, dist is zero)
      beta.mntd.weighted = as.matrix(picante::comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T,exclude.conspecifics=F));
      beta.mntd.weighted
      resultsMatrix[com1,com2]=beta.mntd.weighted[2,1]
    }
    # Adjusting format of results matrix, probably very inefficiently
    resultsMatrixFinal= matrix(as.matrix(resultsMatrix,nrows=length(sampleNames)),ncol=length(sampleNames))
    colnames(resultsMatrixFinal)=sampleNames
    row.names(resultsMatrixFinal)=sampleNames
    resultsMatrixFinal = data.frame(resultsMatrixFinal)
    bMNTD = as.matrix(resultsMatrixFinal)
    row.names(bMNTD) = colnames(bMNTD)
    bMNTD = melt_dist(bMNTD, dist_name = "bMNTD")
    # Add columns with run info
    bMNTD$TimesMixed = i
    bMNTD.vs.vort = rbind(bMNTD.vs.vort,bMNTD)
}

dim(bMNTD.vs.vort) 
head(bMNTD.vs.vort)
# Looks ok...should be 3900 observations (all combos of the 40 samples at each mixing rate (including vortex controls) x 5 treatments)

# Adjust column names 
colnames(bMNTD.vs.vort)
colnames(bMNTD.vs.vort)[1:2] = c("com1","com2")

# Make sure all columns are in the right class
str(bMNTD.vs.vort)
bMNTD.vs.vort$bMNTD = as.numeric(bMNTD.vs.vort$bMNTD)

# IMPORTANT: We will use the same bMNTD.null dataframe and the same BC dataframe (which contains all possible comparisons), created above.
dim(bMNTD.null) #should be 496 total (all combinations for 32 communities in 1x)
head(bMNTD.null)
str(BC)

#### MERGE THE BMNTD AND BC DISSIMILARILY DF'S ####
# Make a column listing the comparison, in each df
head(bMNTD.vs.vort)
bMNTD.vs.vort$Comparison = paste(bMNTD.vs.vort$com1,bMNTD.vs.vort$com2)
bMNTD.null$Comparison = paste(bMNTD.null$com1,bMNTD.null$com2)
BC$Comparison = paste(BC$com1,BC$com2)
            
# Make sure every bMNTD comparison is in the BC dataframe
# (sum should equal first number in dim)
dim(bMNTD.vs.vort)
sum(bMNTD.vs.vort$Comparison %in% BC$Comparison)

dim(bMNTD.null)
sum(bMNTD.null$Comparison %in% BC$Comparison)

# and now, Merge!
bMNTD.BC.vs.Vort = merge(bMNTD.vs.vort, BC, by.x="Comparison", by.y="Comparison", all.x=TRUE, all.y=FALSE)
bMNTD.BC.vs.Vort.null = merge(bMNTD.null, BC, by.x="Comparison", by.y="Comparison", all.x=TRUE, all.y=FALSE)

# Check it out
tail(bMNTD.BC.vs.Vort)
dim(bMNTD.BC.vs.Vort)

tail(bMNTD.BC.vs.Vort.null)
dim(bMNTD.BC.vs.Vort.null)

### Clean up the dataframes
# Get rid of 'artifact' X's, for neatness sake
bMNTD.BC.vs.Vort$com1 = bMNTD.BC.vs.Vort$com1.x
bMNTD.BC.vs.Vort$com2 = bMNTD.BC.vs.Vort$com2.x
bMNTD.BC.vs.Vort$com1 <- gsub("^.{0,1}", "", bMNTD.BC.vs.Vort$com1)
bMNTD.BC.vs.Vort$com2 <- gsub("^.{0,1}", "", bMNTD.BC.vs.Vort$com2)
bMNTD.BC.vs.Vort$Comparison <- paste(bMNTD.BC.vs.Vort$com1, bMNTD.BC.vs.Vort$com2)

# Delete the .x and .y columns
# First, confirm which columns are com1.x, com2.x, com1.y, com2.y. (Do not count the row name column)
head(bMNTD.BC.vs.Vort)
bMNTD.BC.vs.Vort2 = bMNTD.BC.vs.Vort[-c(2, 3, 6,7)]
# confirm that .x and .y columns were dropped, and Bray_Curtis column is still there
head(bMNTD.BC.vs.Vort2)
# All good, reset bMNTD.BC.vs.Vort df
bMNTD.BC.vs.Vort = bMNTD.BC.vs.Vort2

# Do the same for the null dataset
# Get rid of 'artifact' X's, for neatness sake
bMNTD.BC.vs.Vort.null$com1 = bMNTD.BC.vs.Vort.null$com1.x
bMNTD.BC.vs.Vort.null$com2 = bMNTD.BC.vs.Vort.null$com2.x
bMNTD.BC.vs.Vort.null$com1 <- gsub("^.{0,1}", "", bMNTD.BC.vs.Vort.null$com1)
bMNTD.BC.vs.Vort.null$com2 <- gsub("^.{0,1}", "", bMNTD.BC.vs.Vort.null$com2)
bMNTD.BC.vs.Vort.null$Comparison <- paste(bMNTD.BC.vs.Vort.null$com1, bMNTD.BC.vs.Vort.null$com2)

# Delete the .x and .y columns
# First, confirm which columns are com1.x, com2.x, com1.y, com2.y. (Do not count the row name column)
head(bMNTD.BC.vs.Vort.null)
bMNTD.BC.vs.Vort.null2 = bMNTD.BC.vs.Vort.null[-c(2, 3, 6,7,8, 9)]
# confirm that .x and .y columns were dropped, and Bray_Curtis column is still there
head(bMNTD.BC.vs.Vort.null2)
# All good, reset bMNTD.BC.vs.Vort.null df
bMNTD.BC.vs.Vort.null = bMNTD.BC.vs.Vort.null2


# Make sure all columns are in the right class and order
str(bMNTD.BC.vs.Vort)
##if not:
# bMNTD.BC.vs.Vort$bMNTD = as.numeric(bMNTD.BC.vs.Vort$bMNTD)
# bMNTD.BC.vs.Vort$Bray_Curtis = as.numeric(bMNTD.BC.vs.Vort$Bray_Curtis)
# bMNTD.BC.vs.Vort$Rep = as.character(bMNTD.BC.vs.Vort$Rep)
# bMNTD.BC.vs.Vort$TimesMixed = as.character(bMNTD.BC.vs.Vort$TimesMixed)
bMNTD.BC.vs.Vort$TimesMixed = ordered(bMNTD.BC.vs.Vort$TimesMixed,levels=c("1","2","4","8","16","32")) #for graphing later

str(bMNTD.BC.vs.Vort.null)
# bMNTD.BC.vs.Vort.null$bMNTD = as.numeric(bMNTD.BC.vs.Vort.null$bMNTD)
# bMNTD.BC.vs.Vort.null$Bray_Curtis = as.numeric(bMNTD.BC.vs.Vort.null$Bray_Curtis)
# bMNTD.BC.vs.Vort.null$TimesMixed = as.character(bMNTD.BC.vs.Vort.null$TimesMixed)


#### Assign quantiles (for bMNTD and BC) ####

## Read in metadata
meta = read.csv("sample_metadata_connors_fall2018.2021.csv",
                header = TRUE)

head(bMNTD.BC.vs.Vort)
head(meta)
meta = meta[,-c(3,4,5,6,9, 10,11,12,13,14,15)]
colnames(bMNTD.BC.vs.Vort)[colnames(bMNTD.BC.vs.Vort) == "com1"] = "SampleID"


bMNTD.BC.vs.Vort.interm = merge(bMNTD.BC.vs.Vort, meta, by = "SampleID")
head(bMNTD.BC.vs.Vort.interm)
colnames(bMNTD.BC.vs.Vort.interm)= c("com1","Comparison","bMNTD","TimesMixed","Bray_Curtis","SampleID","Rep.com1","VortexControl.com1","GroupOrAgit.com1","Group2.com1")
bMNTD.BC.vs.Vort.interm = merge(bMNTD.BC.vs.Vort.interm, meta, by = "SampleID")
head(bMNTD.BC.vs.Vort.interm)

colnames(bMNTD.BC.vs.Vort.interm)= c("com2","com1","Comparison","bMNTD","TimesMixed","Bray_Curtis",
                             "Rep.com1","VortexControl.com1","GroupOrAgit.com1","Group2.com1",
                             "Rep.com2","VortexControl.com2","GroupOrAgit.com2","Group2.com2")
# Make sorting columns
bMNTD.BC.vs.Vort.interm$Pooled = ifelse(bMNTD.BC.vs.Vort.interm$Group2.com1==bMNTD.BC.vs.Vort.interm$Group2.com2, "Pooled", "Not pooled")
bMNTD.BC.vs.Vort.interm$Vortex = ifelse(bMNTD.BC.vs.Vort.interm$GroupOrAgit.com1=="Agit" & bMNTD.BC.vs.Vort.interm$GroupOrAgit.com2=="Agit", "Both Vortex","Not both vortex")
bMNTD.BC.vs.Vort.interm$Vortex2 = ifelse(bMNTD.BC.vs.Vort.interm$GroupOrAgit.com1=="Agit" | bMNTD.BC.vs.Vort.interm$GroupOrAgit.com2=="Agit", "At least one vortex","No vortex")
head(bMNTD.BC.vs.Vort.interm)
# Make a subset of pooled samples vs vortex contols of the same level of mixing (i.e., Pooled vs. Vortex, or "PvV")
# these are the ones designated as "Not pooled" AND "Not both vortex"
bMNTD.BC.vs.Vort.PvV = subset(bMNTD.BC.vs.Vort.interm, Pooled=="Not pooled" & Vortex =="Not both vortex" & Vortex2=="At least one vortex")

# Make Proc.ID dataframe as a duplicate of the bMNTD.BC.vs.Vort df
Proc.ID = data.frame(bMNTD.BC.vs.Vort.PvV)

# create variable for number of comparisons; this will be used to calculate quantiles
obs = nrow(bMNTD.BC.vs.Vort.null)
# create a vector with even quantile steps
q = seq(0, ((obs-1)/obs), by=(1/obs))
# Assign quantiles (non-parametrically) to each comparison, based on bMNTD values
# first sort rows/observations by bMNTD values
bMNTD.BC.vs.Vort.null = bMNTD.BC.vs.Vort.null[order(bMNTD.BC.vs.Vort.null$bMNTD),]
# then assign quantiles to ascending bMNTD
bMNTD.BC.vs.Vort.null$bMNTD.quant = q
# Similarly, assign quantiles (non-paramentrically) to each comparison, based on Bray_Curtis values
bMNTD.BC.vs.Vort.null = bMNTD.BC.vs.Vort.null[order(bMNTD.BC.vs.Vort.null$Bray_Curtis),]
bMNTD.BC.vs.Vort.null$B_C.quant = q
  
#### Using the null set quantiles, assign quantiles for within-treatment (within mixing set) comparisons
head(bMNTD.BC.vs.Vort.PvV)
# add bMNTD.quant and B_C.quant columns to bMNTD.BC.vs.Vort.PvV dataframe
# which contains all of the within-mixing set comparisons for each treatment
bMNTD.BC.vs.Vort.PvV$bMNTD.quant <- NA
bMNTD.BC.vs.Vort.PvV$B_C.quant <- NA

# Create function that takes a number (from treatment df "bMNTD.BC.vs.Vort.PvV"), and identifies the closest value 
# in the null set df "bMNTD.BC.vs.Vort.null" and populates the trt df (bMNTD.BC.vs.Vort.PvV) with corresponding quantile
closest <- function(reflist, trtvalue){
  reflist[which(abs(reflist-trtvalue)==min(abs(reflist-trtvalue)))]}

# find bMNTD quantiles
for(r in 1:(nrow(bMNTD.BC.vs.Vort.PvV))) {
  bMNTD.BC.vs.Vort.PvV$bMNTD.quant[r] = bMNTD.BC.vs.Vort.null[bMNTD.BC.vs.Vort.null$bMNTD==closest(bMNTD.BC.vs.Vort.null$bMNTD[],bMNTD.BC.vs.Vort.PvV$bMNTD[r]),]$bMNTD.quant
}
# find Bray_Curtis quantiles
for(r in 1:(nrow(bMNTD.BC.vs.Vort.PvV))) {
  bMNTD.BC.vs.Vort.PvV$B_C.quant[r] = bMNTD.BC.vs.Vort.null[bMNTD.BC.vs.Vort.null$Bray_Curtis==closest(bMNTD.BC.vs.Vort.null$Bray_Curtis[],bMNTD.BC.vs.Vort.PvV$Bray_Curtis[r]),]$B_C.quant
}

### Then, assign cutoff quantiles to community assembly processes. 
#create empty column for the process ID
bMNTD.BC.vs.Vort.PvV$Process <- NA
# establish the lower and upper cutoffs, approx equivalent to 0.025th and 0.975th quantiles.
# We could set the cutoffs equal to the target quantile values, but depending on the number of comparisons,
# this won't be a round number; and, we have to account for the 1st quantile step = 0!
# so we will use "<=" for lower quantile, to capture the 0th quantile
# and ">" for upper quantile (not >=).

# Find the closest element of q to 0.05, i.e. how many elements get us closest to 0.05
howmany= which(abs(q-0.05)==min(abs(q-0.05)))
# Now make it an even number
howmany2 = if((howmany %% 2) == 0) {
  which(abs(q-0.05)==min(abs(q-0.05)))
} else {
  (which(abs(q-0.05)==min(abs(q-0.05))) - 1) # Could add 1 instead..
}
# Establish the lower quantile cutoff
lowerquant= q[howmany2/2]
# establish the upper quantile cutoff.
upperquant= q[nobs(q)-(howmany2/2)]

# assign process ID's, considering both bMNTD and BC-Dis cutoffs for each comparison
  for(r in 1:(nrow(bMNTD.BC.vs.Vort.PvV))) {
    if (bMNTD.BC.vs.Vort.PvV$bMNTD.quant[r] <= lowerquant) {
      bMNTD.BC.vs.Vort.PvV$Process[r] = "Homogeneous selection"
    } else if (bMNTD.BC.vs.Vort.PvV$B_C.quant[r] <= lowerquant && bMNTD.BC.vs.Vort.PvV$bMNTD.quant[r] > lowerquant && bMNTD.BC.vs.Vort.PvV$bMNTD.quant[r] <= upperquant) {
      bMNTD.BC.vs.Vort.PvV$Process[r] = "Homogenizing dispersal"
      } else if (bMNTD.BC.vs.Vort.PvV$bMNTD.quant[r] > upperquant) {
      bMNTD.BC.vs.Vort.PvV$Process[r] = "Variable selection"
    } else if (bMNTD.BC.vs.Vort.PvV$B_C.quant[r] > upperquant && bMNTD.BC.vs.Vort.PvV$bMNTD.quant[r] > lowerquant && bMNTD.BC.vs.Vort.PvV$bMNTD.quant[r] <= upperquant) {
      bMNTD.BC.vs.Vort.PvV$Process[r] = "Dispersal limitation"
    } else if (bMNTD.BC.vs.Vort.PvV$B_C.quant[r] > lowerquant && bMNTD.BC.vs.Vort.PvV$B_C.quant[r] <= upperquant && bMNTD.BC.vs.Vort.PvV$bMNTD.quant[r] > lowerquant && bMNTD.BC.vs.Vort.PvV$bMNTD.quant[r] <= upperquant) {
      bMNTD.BC.vs.Vort.PvV$Process[r] = "Undominated"
    } else bMNTD.BC.vs.Vort.PvV$Process[r] = NA
  }
  
head(bMNTD.BC.vs.Vort.PvV)


#### TALLY/SUMMARIZE THE COMMUNITY ASSEMBLY PROCESSES ASSIGNED TO WITHIN-TREATMENT COMPARISSONS ####

# Check that no ID assignments were missed; there should be no NAs for any of the Iter. columns
bMNTD.BC.vs.Vort.PvV %>%
  select_if(function(x) any(is.na(x))) %>%
  summarise_each(funs(sum(is.na(.))))

### Reorganize data for summarizing/formatting for stats and graphing

# Summarize the process ID assignments.
ProcessSummary = bMNTD.BC.vs.Vort.PvV %>%
  group_by(TimesMixed,Process)%>%
  dplyr::summarize(Count = n())
str(ProcessSummary)


#### STACKED BAR GRAPH--TREATMENT COMMUNITIES ####
# Rename dataframe, for ease of use
ProcID2 = data.frame(ProcessSummary)
head(ProcID2)
str(ProcID2)
ProcID2
#
# Order factors for graphing
ProcID2$TimesMixed = ordered(ProcID2$TimesMixed,levels=c("2", "4", "8", "16", "32"))
ProcID2$Process = ordered(ProcID2$Process, levels=c("Undominated",
                                            "Homogenizing dispersal",
                                            "Homogeneous selection",
                                            "Dispersal limitation",
                                            "Variable selection"))


### Stacked bar graph with error bars
ID_Palette = c("#CCCCCC", "#336699", "#99CC99", "#FFCC66", "#CC3333")
n.obs = 256
plot1 = ggplot(ProcID2, aes(x=factor(TimesMixed), y=Count/n.obs)) +
  geom_bar(position = "stack", stat="identity",aes(fill=Process)) +
  expand_limits(y=c(0,1.04)) +
  scale_fill_manual(values = ID_Palette)+
  scale_x_discrete(labels=c("2"="2\u00d7", "4"="4\u00d7", "8"="8\u00d7", "16"="16\u00d7", "32"="32\u00d7")) +
  labs(x = "Mixing frequency \nComparisons between pooled tubes + vortex controls", y = "Proportion of comparisons", fill = "Community assembly process") +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+
  #theme(legend.position = "none") 
plot1 #600 x 500
