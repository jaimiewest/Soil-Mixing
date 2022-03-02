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


#### Create randomized sets of the 1X controls  ####
#### and assign quantiles (for bMNTD and BC) based on comparisons within those sets

## Read in metadata
meta = read.csv("sample_metadata_connors_fall2018.2021.csv",
                header = TRUE)

# Make a subset of metadata, 1x  only
meta.1x = subset(meta, IncubBlank=="N" & ExtrBlank =="N" & TimesMixed =="1")
# Drop unnecessary columns; keep SampleID, TimesMixed, and GroupOrAgit
head(meta.1x)
meta.1x = meta.1x[,c(1, 3, 8)]

### For loop that randomly assigns the 32 1x controls to 4 sets

# Make Process.ID dataframe as a duplicate of the bMNTD.BC df
Process.ID = data.frame(bMNTD.BC)

# # Option to create a dataframe to track the 1x set assignment
# null = data.frame()
# null = data.frame(iter=integer(0), set=integer(0), Comparison=character(0), bMNTD=numeric(0), Bray_Curtis=numeric(0), 
#                   bMNTD.quant=numeric(0), B_C.quant=numeric(0))

iterations = 999 # took about 3 mins/999 iterations, and almost 60 mins/2600 iterations
for (i in 1:iterations){
  re = sample(factor(rep(1:4, length.out=nrow(meta.1x)), labels=paste0(1:4)))
  meta.1x$set = re
  # create intermediate df for this iteration
  inter.df = data.frame(iter=integer(0), set=integer(0), Comparison=character(0), bMNTD=numeric(0), Bray_Curtis=numeric(0), 
                        bMNTD.quant=numeric(0), B_C.quant=numeric(0))
  for (j in c(1,2,3,4)){
    # create df where com1 column of bMNTD.BC.null matches rep (1:4) of null set.j
    df.set = merge(subset(meta.1x, set == j, select = c("SampleID", "TimesMixed", "GroupOrAgit", "set")), 
                    bMNTD.BC.null[, c("Comparison", "bMNTD", "Bray_Curtis", "com1", "com2")], 
                    by.x="SampleID", by.y="com1")
    # Further wittle down such that com2 column of bMNTD.BC.null ALSO matches rep (1:4) of null set.j
    df.set = merge(subset(meta.1x, set == j, select = c("SampleID", "TimesMixed", "GroupOrAgit", "set")), 
                    df.set[, c("Comparison", "bMNTD", "Bray_Curtis", "com2")],
                    by.x="SampleID", by.y="com2")
    # Add in each rep to intermediate df
    inter.df = Reduce(
      function(...) merge(..., all=TRUE), list(inter.df, df.set))
  }
  #record the iteration number
  inter.df$iter = i
  
  #### Assign quantiles for the null sets 
  
  # create variable for number of comparisons; this will be used to calculate quantiles
  obs = nrow(inter.df)
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
  # (thus allowing for dual assignment)
  for(r in 1:(nrow(bMNTD.BC))) {
    if (bMNTD.BC$B_C.quant[r] <= lowerquant && bMNTD.BC$bMNTD.quant[r] <= lowerquant) {
      bMNTD.BC$Process[r] = "Homogeneous Selection & Homogenizing Dispersal"
    } else if (bMNTD.BC$B_C.quant[r] <= lowerquant && bMNTD.BC$bMNTD.quant[r] > lowerquant && bMNTD.BC$bMNTD.quant[r] <= upperquant) {
      bMNTD.BC$Process[r] = "Homogenizing Dispersal"
    } else if (bMNTD.BC$B_C.quant[r] > lowerquant && bMNTD.BC$B_C.quant[r] <= upperquant && bMNTD.BC$bMNTD.quant[r] <= lowerquant) {
      bMNTD.BC$Process[r] = "Homogeneous Selection"
    } else if (bMNTD.BC$B_C.quant[r] > upperquant && bMNTD.BC$bMNTD.quant[r] > upperquant) {
      bMNTD.BC$Process[r] = "Dispersal Limitation & Variable Selection"
    } else if (bMNTD.BC$B_C.quant[r] > upperquant && bMNTD.BC$bMNTD.quant[r] > lowerquant && bMNTD.BC$bMNTD.quant[r] <= upperquant) {
      bMNTD.BC$Process[r] = "Dispersal Limitation"
    } else if (bMNTD.BC$B_C.quant[r] > lowerquant && bMNTD.BC$B_C.quant[r] <= upperquant && bMNTD.BC$bMNTD.quant[r] > upperquant) {
      bMNTD.BC$Process[r] = "Variable Selection"
    } else if (bMNTD.BC$B_C.quant[r] > lowerquant && bMNTD.BC$B_C.quant[r] <= upperquant && bMNTD.BC$bMNTD.quant[r] > lowerquant && bMNTD.BC$bMNTD.quant[r] <= upperquant) {
      bMNTD.BC$Process[r] = "Undominated"
    } else bMNTD.BC$Process[r] = NA
  }
  
  ### Add this iteration's process ID column to the Process.ID dataframe.
  Process.ID$Process = bMNTD.BC$Process
  ### Give the Process column a unique name for this iteration
  #names(bMNTD.BC)[ncol(bMNTD.BC)]= paste0("Process.Iter.",i)
  names(bMNTD.BC)[ncol(bMNTD.BC)]= paste0("Iter.",i,sep="")
  
  # ## OPTIONAL: Add this iteration, with assigned quantiles, to the null df
  # ## ...in case you want to check the process; no need to record or save this, though.
  # null = Reduce(
  # function(...) merge(..., all=TRUE), list(null, inter.df[-5]))
}

## Option to save data
write.csv(bMNTD.BC,"Derived_data/ProcessID_Assignments_IterIs999.csv")


#### TALLY/SUMMARIZE THE COMMUNITY ASSEMBLY PROCESSES ASSIGNED TO WITHIN-TREATMENT COMPARISONS ####

# Check that no ID assignments were missed; there should be no NAs for any of the Iter. columns
bMNTD.BC %>%
  select_if(function(x) any(is.na(x))) %>%
  summarise_each(funs(sum(is.na(.))))

### Reorganize data for summarizing/formatting for stats and graphing

# Expand df to include a column that holds the iteration number and then just one column 
# where all the processes go; summarize from there.

# Rearrange the data.
ProcessGathered = gather(data=bMNTD.BC, key=Iteration, value = Process, Iter.1:(paste(colnames(bMNTD.BC)[length(bMNTD.BC)])),factor_key=TRUE)
# head(ProcessGathered)
# tail(ProcessGathered)

# Double check to be sure there are no NAs under Process (due to miscoding of process IDs); Process should return FALSE
apply(ProcessGathered, 2, function(x) any(is.na(x)))

# Check to be sure the number of rows in ProcessGathered = number of rows in bMNTD*iterations; i.e. TRUE
dim(ProcessGathered)[1] == (dim(bMNTD.BC)[1]) * iterations

# Now, we can summarize the process ID assignments.
ProcessSummary = ProcessGathered %>%
  group_by(TimesMixed,VortexControl,Iteration,Process)%>%
  dplyr::summarize(Count = n())
str(ProcessSummary)

# Use spread, the opposite of gather, to reformat data back into 'wide' format:
ProcessSummarySpread = spread(data=ProcessSummary,key = Process,value = Count, fill=0)
head(ProcessSummarySpread)

# Option to save data
write.csv(ProcessSummarySpread, "Derived_data/Community_Assembly_Process_ID_Summary.csv")

#### STACKED BAR GRAPH--TREATMENT COMMUNITIES ####

# Rename dataframe, for ease of use
ProcID = data.frame(ProcessSummarySpread)

# Add column with number of observation for each row
ProcID$n.obs=rowSums(ProcID[,4:10])

### Subset for TREATMENT communities
ProcID.trt = subset(ProcID, VortexControl=="N")
ProcID.trt

# Find the mean proportion of each process ID, within treatments
ProcID.prop = ProcID.trt %>% 
  group_by(TimesMixed) %>%
  dplyr::summarise(n=n(),
            "Variable Selection" = mean(Variable.Selection/n.obs),
            "Homogeneous Selection" = mean(Homogeneous.Selection/n.obs),
            "Homogenizing Dispersal" = mean(Homogenizing.Dispersal/n.obs),
            "Homogeneous Selection & Homogenizing Dispersal" = mean(Homogeneous.Selection...Homogenizing.Dispersal/n.obs),
            "Dispersal Limitation & Variable Selection" = mean(Dispersal.Limitation...Variable.Selection/n.obs),
            "Dispersal Limitation" = mean(Dispersal.Limitation/n.obs),
            "Undominated" = mean(Undominated/n.obs)) %>%
  gather("ID", "proportion", - c(TimesMixed, n), factor_key=TRUE)  
head(ProcID.prop)


# Find the standard deviation of the proportions
ProcID.SD = ProcID.trt %>% 
  group_by(TimesMixed) %>%
  dplyr::summarise(n=n(),
                   "Variable Selection" = sd(Variable.Selection/n.obs),
                   "Homogeneous Selection" = sd(Homogeneous.Selection/n.obs),
                   "Homogenizing Dispersal" = sd(Homogenizing.Dispersal/n.obs),
                   "Homogeneous Selection & Homogenizing Dispersal" = sd(Homogeneous.Selection...Homogenizing.Dispersal/n.obs),
                   "Dispersal Limitation & Variable Selection" = sd(Dispersal.Limitation...Variable.Selection/n.obs),
                   "Dispersal Limitation" = sd(Dispersal.Limitation/n.obs),
                   "Undominated" = sd(Undominated/n.obs)) %>%
  gather("ID", "SD", - c(TimesMixed, n), factor_key=TRUE)
head(ProcID.SD)

# Join proportion and SD data together
ProcID.2 = merge(ProcID.prop,ProcID.SD, by=c("TimesMixed","n","ID"))

# Order factors for graphing
ProcID.2$TimesMixed = ordered(ProcID.2$TimesMixed,levels=c("2", "4", "8", "16", "32"))
ProcID.2$ID = ordered(ProcID.2$ID, levels=c("Undominated",
                                            "Homogenizing Dispersal",
                                            "Homogeneous Selection & Homogenizing Dispersal",
                                            "Homogeneous Selection",
                                            "Variable Selection",
                                            "Dispersal Limitation & Variable Selection",
                                            "Dispersal Limitation"))

# Establish y position so that our error bars are positioned relative to stacking
ProcID.2$y_pos = NA
ProcID.2$y_pos[ProcID.2$ID == "Dispersal Limitation"] =
  ProcID.2$proportion[ProcID.2$ID =="Dispersal Limitation"]
ProcID.2$y_pos[ProcID.2$ID == "Dispersal Limitation & Variable Selection"] =
  ProcID.2$proportion[ProcID.2$ID =="Dispersal Limitation & Variable Selection"] + ProcID.2$y_pos[ProcID.2$ID =="Dispersal Limitation"]
ProcID.2$y_pos[ProcID.2$ID == "Variable Selection"] =
  ProcID.2$proportion[ProcID.2$ID =="Variable Selection"] + ProcID.2$y_pos[ProcID.2$ID =="Dispersal Limitation & Variable Selection"]
ProcID.2$y_pos[ProcID.2$ID == "Homogeneous Selection"] =
  ProcID.2$proportion[ProcID.2$ID =="Homogeneous Selection"] + ProcID.2$y_pos[ProcID.2$ID =="Variable Selection"]
ProcID.2$y_pos[ProcID.2$ID == "Homogeneous Selection & Homogenizing Dispersal"] =
  ProcID.2$proportion[ProcID.2$ID =="Homogeneous Selection & Homogenizing Dispersal"] + ProcID.2$y_pos[ProcID.2$ID =="Homogeneous Selection"]
ProcID.2$y_pos[ProcID.2$ID == "Homogenizing Dispersal"] =
  ProcID.2$proportion[ProcID.2$ID =="Homogenizing Dispersal"] + ProcID.2$y_pos[ProcID.2$ID =="Homogeneous Selection & Homogenizing Dispersal"]
ProcID.2$y_pos[ProcID.2$ID == "Undominated"] =
  ProcID.2$proportion[ProcID.2$ID =="Undominated"] + ProcID.2$y_pos[ProcID.2$ID =="Homogenizing Dispersal"]
head(ProcID.2)

### Stacked bar graph with error bars
dodge = position_dodge2(width = 0.3, padding = 0.3)
ID_Palette = c("#CCCCCC", "#0072B2", "#009999", "#006633", "#FFCC00", "#FF9900", "#CC3300")
Trt.plot = ggplot(ProcID.2[which(ProcID.2$proportion>0),], aes(x=factor(TimesMixed), y=proportion)) +
  geom_bar(position = "stack", stat="identity",aes(fill=ID)) +
  scale_fill_manual(values = ID_Palette) +
  geom_errorbar(aes(ymin= y_pos - SD, ymax= y_pos + SD), width=0.2, position = dodge)+#position_dodge(0.1)) + # 1SD errr bars
  scale_x_discrete(labels=c("2"="2\u00d7", "4"="4\u00d7", "8"="8\u00d7", "16"="16\u00d7", "32"="32\u00d7")) +
  labs(x = "Mixing frequency", y = "Proportion of comparisons", fill = "Community assembly process") +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  theme_bw()
Trt.plot

#### STACKED BAR GRAPH--VORTEX CONTROL COMMUNITIES ####

### Subset for Vortex CONTROL communities
ProcID.Vort = subset(ProcID, VortexControl=="Y")

# Find the mean proportion of each process ID, within TimesMixed
ProcID.Vort.prop = ProcID.Vort %>% 
  group_by(TimesMixed) %>%
  dplyr::summarise(n=n(),
            "Variable Selection" = mean(Variable.Selection/n.obs),
            "Homogeneous Selection" = mean(Homogeneous.Selection/n.obs),
            "Homogenizing Dispersal" = mean(Homogenizing.Dispersal/n.obs),
            "Homogeneous Selection & Homogenizing Dispersal" = mean(Homogeneous.Selection...Homogenizing.Dispersal/n.obs),
            "Dispersal Limitation & Variable Selection" = mean(Dispersal.Limitation...Variable.Selection/n.obs),
            "Dispersal Limitation" = mean(Dispersal.Limitation/n.obs),
            "Undominated" = mean(Undominated/n.obs)) %>%
  gather("ID", "proportion", - c(TimesMixed, n), factor_key=TRUE)  
str(ProcID.Vort.prop)

# Find the standard deviation of the proportions
ProcID.Vort.SD = ProcID.Vort %>%
  group_by(TimesMixed) %>%
  dplyr::summarise(n=n(),
                   "Variable Selection" = sd(Variable.Selection/n.obs),
                   "Homogeneous Selection" = sd(Homogeneous.Selection/n.obs),
                   "Homogenizing Dispersal" = sd(Homogenizing.Dispersal/n.obs),
                   "Homogeneous Selection & Homogenizing Dispersal" = sd(Homogeneous.Selection...Homogenizing.Dispersal/n.obs),
                   "Dispersal Limitation & Variable Selection" = sd(Dispersal.Limitation...Variable.Selection/n.obs),
                   "Dispersal Limitation" = sd(Dispersal.Limitation/n.obs),
                   "Undominated" = sd(Undominated/n.obs)) %>%
  gather("ID", "SD", - c(TimesMixed, n), factor_key=TRUE)
head(ProcID.Vort.SD)

# Join proportion and SD data together
ProcID.Vort.2 = merge(ProcID.Vort.prop,ProcID.Vort.SD, by=c("TimesMixed","n","ID"))

# Order factors for graphing
ProcID.Vort.2$TimesMixed = ordered(ProcID.Vort.2$TimesMixed,levels=c("2", "4", "8", "16", "32"))
ProcID.Vort.2$ID = ordered(ProcID.Vort.2$ID, levels=c("Undominated",
                                                        "Homogenizing Dispersal",
                                                        "Homogeneous Selection & Homogenizing Dispersal",
                                                        "Homogeneous Selection",
                                                        "Variable Selection",
                                                        "Dispersal Limitation & Variable Selection",
                                                        "Dispersal Limitation"))


# Establish y position so that our error bars are positioned relative to stacking
ProcID.Vort.2$y_pos = NA
ProcID.Vort.2$y_pos[ProcID.Vort.2$ID == "Dispersal Limitation"] =
  ProcID.Vort.2$proportion[ProcID.Vort.2$ID =="Dispersal Limitation"]
ProcID.Vort.2$y_pos[ProcID.Vort.2$ID == "Dispersal Limitation & Variable Selection"] =
  ProcID.Vort.2$proportion[ProcID.Vort.2$ID =="Dispersal Limitation & Variable Selection"] + ProcID.Vort.2$y_pos[ProcID.Vort.2$ID =="Dispersal Limitation"]
ProcID.Vort.2$y_pos[ProcID.Vort.2$ID == "Variable Selection"] =
  ProcID.Vort.2$proportion[ProcID.Vort.2$ID =="Variable Selection"] + ProcID.Vort.2$y_pos[ProcID.Vort.2$ID =="Dispersal Limitation & Variable Selection"]
ProcID.Vort.2$y_pos[ProcID.Vort.2$ID == "Homogeneous Selection"] =
  ProcID.Vort.2$proportion[ProcID.Vort.2$ID =="Homogeneous Selection"] + ProcID.Vort.2$y_pos[ProcID.Vort.2$ID =="Variable Selection"]
ProcID.Vort.2$y_pos[ProcID.Vort.2$ID == "Homogeneous Selection & Homogenizing Dispersal"] =
  ProcID.Vort.2$proportion[ProcID.Vort.2$ID =="Homogeneous Selection & Homogenizing Dispersal"] + ProcID.Vort.2$y_pos[ProcID.Vort.2$ID =="Homogeneous Selection"]
ProcID.Vort.2$y_pos[ProcID.Vort.2$ID == "Homogenizing Dispersal"] =
  ProcID.Vort.2$proportion[ProcID.Vort.2$ID =="Homogenizing Dispersal"] + ProcID.Vort.2$y_pos[ProcID.Vort.2$ID =="Homogeneous Selection & Homogenizing Dispersal"]
ProcID.Vort.2$y_pos[ProcID.Vort.2$ID == "Undominated"] =
  ProcID.Vort.2$proportion[ProcID.Vort.2$ID =="Undominated"] + ProcID.Vort.2$y_pos[ProcID.Vort.2$ID =="Homogenizing Dispersal"]


### Graph: stacked bar with error bars
dodge = position_dodge2(width = 0.3, padding = 0.3)
#ID_Palette = c("#CCCCCC", "#0072B2", "#009999", "#006633", "#FFCC00", "#FF9900", "#CC3300") #full set for all process IDs
ID_Palette = c("#CCCCCC", "#0072B2", "#006633", "#FFCC00", "#FF9900", "#CC3300")
Vort.plot = ggplot(ProcID.Vort.2[which(ProcID.Vort.2$proportion>0),], aes(x=factor(TimesMixed), y=proportion, fill=ID)) +
  geom_bar(position = "stack", stat="identity") +
  scale_fill_manual(values = ID_Palette) +
  geom_errorbar(aes(ymin= y_pos - SD, ymax= y_pos + SD), width=0.2, position = dodge) + # SD error bars
  scale_x_discrete(labels=c("2"="2\u00d7", "4"="4\u00d7", "8"="8\u00d7", "16"="16\u00d7", "32"="32\u00d7")) +
  labs(x = "Mixing frequency, vortex controls", y = "Proportion of comparisons", fill = "Community assembly process") +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  scale_y_continuous(breaks=seq(0,1,0.2)) +
  theme_bw()
Vort.plot

