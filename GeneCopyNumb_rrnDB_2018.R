library(Biostrings)
library(phyloseq)
library(plyr)
library(dplyr)
library(vegan)
library(ggplot2)
library(DescTools) #for Dunnett's test

# We want to test how the mean predicted 16S copy number changes with frequency of soil mixing.

#### Adding RDP Classifier data ####

# From the rrnDB website: "Estimate is an on-line interface to the RDP Classifier tool, including adjustment of 
# relative abundance of taxons based on 16S gene copy number data from rrnDB."
# "Estimate runs RDP Classifier version 2.12 using 16S training set #16 incorporating current rrnDB copy number data. 
# The necessary training files were re-created following the documented use of RDP Classifier's 'train' command, 
# replacing the copy number file from the original training set with one derived from the most recent downloadable 
# pan-taxa statistics."

# 1. https://rrndb.umms.med.umich.edu/estimate/run_classifier
# 2. Navigate to the “Download” page
# 3. Download “rrnDB-5.7_pantaxa_stats_RDP.tsv”, or the most recent version of this
# 4. Navigate back to rrnDB "Estimate" page and upload the fasta file for the full dataset (confidence cutoff 0.8)
# 5. Wait several minutes for it to run
# 6. Download the Classification assignment file (“dna-sequences.tsv”)


# Load in sequence classifications file
RDP = read.csv("rrnDB/CON18_dna-sequences.tsv",header=FALSE,sep=";")
head(RDP)
# We want to extract the genus they assigned for each OTU.
# Create function to extract genus, if present
GenusGenerator <- function(taxon) {
  Genus = strsplit(gsub(".*family\\s*|genus*", "", taxon),split="\t")[[1]][2]
  return(Genus)
}

# Extract the genus
RDP$GenusRDP = sapply(RDP$V1,GenusGenerator)
head(RDP)

# Might as well pull out OTU ID to be certain
OTUGenerator <- function(taxon) {
  OTU = strsplit(paste(taxon),split="\t")[[1]][1]
  return(OTU)
}
RDP$OTU = sapply(RDP$V1,OTUGenerator)
head(RDP)

# Ok, we've got what we need.
# Trim it down
RDP = RDP[,c(2,3)]

# Can now pull data from RRNDB.
# Reading in the rrnDB v5.5 file
rrnDB = read.csv("rrnDB/rrnDB-5.7_pantaxa_stats_RDP.tsv",sep="\t")
head(rrnDB)

# Creating a list of genera in the DB
rrnDBGenera = as.character(rrnDB[rrnDB$rank=="genus",]$name)


# Matching up genus name with mean predicted copy number
for (i in 1:length(RDP$GenusRDP)){
  GenusRDP = paste(RDP$GenusRDP[i])
  CopyNum = ifelse(GenusRDP %in% rrnDBGenera, rrnDB[rrnDB$name==GenusRDP,9],"")
  RDP$CopyNum[i] = CopyNum
}
tail(RDP)

### OPTIONAL--test for removal of Nocardioides gene copy number
#RDP = subset(RDP, GenusRDP!="Nocardioides")

# Bring in ps object and normalize to relative abundances
ps.full <- readRDS("ps.CON18.remake")
ps.norm <- transform_sample_counts(ps.full, function(x) x / sum(x))

# Working with melted phyloseq object
mdf = psmelt(ps.norm)

# Add the rrnDB copy number data to the melted phyloseq object
mdf = plyr::join(mdf,RDP,by="OTU")
mdf$CopyNum = as.numeric(mdf$CopyNum)
mdf$Abundance = as.numeric(mdf$Abundance)

head(mdf)

# From Nemergut et al. (2016) - OTU data were then normalized (standardized) for copy number 
# by dividing by copy number. For each sample, we calculated the community aggregated trait value 
# (weighted mean) by taking the product of the estimated operon copy number and the relative abundance 
# for each OTU, and summing this value across all OTUs in a sample. 

# So, first, we divide abundance by copy number
# Then, we re-calculate the relative abundanace, now adjusted for copy number
# The risk there, is, for any organisms without assigned copy numbers, they are excluded from this calculation.
# To check:

d = mdf %>%
  dplyr::group_by(VortexControl, TimesMixed, Sample)%>%
  dplyr::filter(is.na(CopyNum))%>%
  dplyr::summarize(NoCopyNum = sum(Abundance))
hist(d$NoCopyNum)
d
# Not too bad, but a wide range

## Check out distribution of NoCopyNum across treatments.
# First, mixed soil trts
# create df for just the mixing experiment
# and then another df that excludes the TimesMixed column (for the grey background full dataset)
dat.trt=subset(d, VortexControl=="N" & TimesMixed!="Initial")
head(dat.trt)
dat.notrt=dat.trt[, -2]
p = ggplot(dat.trt, aes(x=NoCopyNum, fill=TimesMixed)) +
  geom_histogram(data = dat.notrt, fill = "grey", alpha = 0.5) +
  geom_histogram(color = "black") +
  facet_wrap(~TimesMixed) +
  theme_bw()
p = p + scale_color_manual(values=palette)
p = p + xlab("Proportion of relative abundance for which OTUs lack a predicted copy number")
p = p + ylab("Count, number of samples")
p 

# There does appear to be a connection between mixing frequency, and copy number data in the rrnDB.
# Likely that some of the high-abundance OTUs in frequently mixed trts are in the database, weighing heavily on this proportion.


# Next, check out Vortex Control trts
# create df for just the vortex controls
# and then another df that excludes the TimesMixed column (for the grey background full dataset)
dat.vort=subset(d, VortexControl=="Y" & TimesMixed!="Initial")
head(dat.trt)
dat.notrt=dat.vort[, -2]
p = ggplot(dat.vort, aes(x=NoCopyNum, fill=TimesMixed)) +
  geom_histogram(data = dat.notrt, fill = "grey", alpha = 0.5) +
  geom_histogram(color = "black") +
  facet_wrap(~TimesMixed) +
  theme_bw()
p = p + scale_color_manual(values=palette)
p = p + xlab("Proportion of relative abundance for which OTUs lack a predicted copy number")
p = p + ylab("Count, number of samples")
p 

# Similar pattern as the mixed soil treatments.


### Optional mean replacement: unassigned taxa are given the mean copy number across dataset
### We did not use this option.
# meanCopyNum = RDP%>%
#   group_by(GenusRDP,CopyNum)%>%
#   dplyr::summarize(N=n())
# meanCopyNum = mean(as.numeric(meanCopyNum$CopyNum),na.rm=TRUE)
# meanCopyNum
# head(mdf)

# Calculating weighted mean copy numbers (with optional mean replacement; not used):
df = mdf %>%
  #dplyr::mutate(CopyNum = ifelse(is.na(CopyNum) | CopyNum =="NA",meanCopyNum,CopyNum))%>%
  dplyr::filter(!is.na(CopyNum))%>%
  dplyr::mutate(WtAbund = Abundance/CopyNum)%>%
  dplyr::group_by(Sample, Rep, TimesMixed, Initial, VortexControl, GroupOrAgit, Group, Group2)%>%
  dplyr::mutate(AdjAbund = WtAbund/sum(WtAbund))%>%
  dplyr::mutate(WtCopyNum = CopyNum*AdjAbund)%>%
  dplyr::summarize(WtMeanCopyNum = sum(WtCopyNum,na.rm=TRUE))
hist(df$WtMeanCopyNum)
head(df)

#write.csv(df, 'Derived_data/Weighted_Mean_16S_Gene_Copy_number.csv')
#write.csv(df, 'Derived_data/Weighted_Mean_16S_Gene_Copy_number_NO_Nocardioides.csv')

#df <- read.csv(file = 'Derived_data/Weighted_Mean_16S_Gene_Copy_number.csv')
#df <- read.csv(file = 'Derived_data/Weighted_Mean_16S_Gene_Copy_number_NO_Nocardioides.csv')


# Plot the results, roughly
p = ggplot(df,aes(x=TimesMixed,y=WtMeanCopyNum,
                     color=VortexControl)) +
  theme_bw() +
  geom_point(size=3,alpha=0.8) +
  labs(x = "Times mixed", 
       y = "Weighted mean predicted\n16S rRNA gene copy number",
       color = "Vortex Control")
p

### Plot the results, neatly
df$TimesMixed = ordered(df$TimesMixed, levels=c("Initial","1","2","4","8","16","32"))
# Group2 is a variable used to tease out each mixing set within mixing treatments.
df$Group2 = ordered(df$Group2, levels=c("Initial","1X","U","V","W","X","Q","R","S","T","M","N","O","P","I","J","K","L","E","F","G","H","ZAgit2","ZAgit4","ZAgit8","ZAgit16","ZAgit32"))

p.2 = ggplot(df,aes(x=TimesMixed,y=WtMeanCopyNum,color=Group2))
p.2 = p.2 + geom_jitter(alpha = 0.1)
p.2 = p.2 + geom_boxplot(alpha = 0.7) + #expand_limits(y=c(0.0,3.0)) +
  labs(x="Mixing frequency", y="Weighted mean predicted\n16S rRNA gene copy number", title=NULL)
p.2 = p.2 + scale_x_discrete(labels=c("Initial"="Initial","1"="1\u00d7 (control)", "2"="2\u00d7", "4"="4\u00d7", "8"="8\u00d7", "16"="16\u00d7", "32"="32\u00d7"))
p.2 = p.2 + scale_color_manual(breaks = c("Initial", "1X","U","Q","M","I","E","ZAgit2"),
                               values = c("gray55", #initial
                                          "#440154FF", #1X
                                          "#404788FF","#404788FF","#404788FF","#404788FF",#2X
                                          "#2D708EFF","#2D708EFF","#2D708EFF","#2D708EFF",#4X 
                                          "#20A387FF","#20A387FF","#20A387FF","#20A387FF",#8X 
                                          "#73D055FF","#73D055FF","#73D055FF","#73D055FF",#16X 
                                          "#FDE725FF","#FDE725FF","#FDE725FF","#FDE725FF",#32X 
                                          "grey20","grey20","grey20","grey20","grey20"),#VortexControls
                               name="", 
                               labels=c("Initial","1\u00d7","2\u00d7","4\u00d7","8\u00d7","16\u00d7","32\u00d7","Vortex Control"))
p.2 = p.2 + theme_bw()
p.2


### Summarize the results

# Calculating means across treatments:
head(df)
df.mean = df %>%
  dplyr::filter(!is.na(CopyNum))%>%
  dplyr::group_by(TimesMixed, VortexControl)%>%
  dplyr::summarize(mean_WtMeanCopyNum=mean(WtMeanCopyNum))
df.mean

head(dat.trt)
CN.mean = dat.trt %>%
  dplyr::group_by(TimesMixed, VortexControl)%>%
  dplyr::summarize(mean_NoCopyNum=mean(NoCopyNum))
CN.mean 



# ANOVA
data.trt = subset(df, VortexControl=="N" & TimesMixed!="Initial" )
head(data.trt)
ano = aov(WtMeanCopyNum ~ TimesMixed, data = data.trt)
summary(ano)
par(mfrow=c(2,2))
plot(ano)
par(mfrow=c(1,1))

# tukey = TukeyHSD(ano, conf.level = 0.95)
# tukey

Dunnett = DunnettTest(data.trt$WtMeanCopyNum, data.trt$TimesMixed)
Dunnett


data.vort = subset(df, VortexControl=="Y" | TimesMixed=="1")
anoVort = aov(WtMeanCopyNum ~ TimesMixed, data = data.vort)
summary(anoVort)
par(mfrow=c(2,2))
plot(anoVort)
par(mfrow=c(1,1))

Dunnett = DunnettTest(data.vort$WtMeanCopyNum, data.vort$TimesMixed)
Dunnett
