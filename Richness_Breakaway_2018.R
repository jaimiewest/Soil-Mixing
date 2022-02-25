# Load required packages

# install.packages("devtools")
# devtools::install_github("adw96/breakaway")
library(breakaway)
# You may need to update Breakaway and then close and reopen R

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
library(phyloseq)
library(tidyverse)
library(plyr)
library(dplyr) # load plyr, then dplyr in that order

# Read in phyloseq object
ps = readRDS("ps.CON18.remake")
ps

# Drop the OTUs that have zeros across all samples
ps.nozero = subset_taxa(ps,taxa_sums(ps)>0)
ps.nozero

# Collect sample data
SamDat = data.frame(sample_data(ps.nozero))

#### Generate Breakaway Richness Estimates ####
#tutorial: https://adw96.github.io/breakaway/articles/breakaway.html

# Set OTU table
otu_data = t(otu_table(ps.nozero))
# Set sample data
meta_data = sample_data(ps.nozero)
# Had to flip OTU table so rownames match sample data
head(colnames(otu_data) == rownames(meta_data))

# Run Breakaway's frequency table list function
frequencytablelist = build_frequency_count_tables(otu_data)

# Check out one of them (#63)
head(frequencytablelist[[63]])

# Try Breakaway on a couple of samples
breakaway(frequencytablelist[[1]])
breakaway(frequencytablelist[[60]])

# Because no plot pops up, we know that we're dealing with the WLRM
# That's because dada2 won't allow singletons to pass through

# Run the richness estimator (breakaway) on all our samples (lists of frequency tables)
RichEsts = lapply(frequencytablelist,breakaway)

# Pull out the estimates, errors, and the model
Estimate = as.matrix(map(RichEsts, "estimate"))
Error = as.matrix(map(RichEsts, "error"))
Model = as.matrix(map(RichEsts, "model"))
df = data.frame(Estimate,Error,Model)

# Add sample ID column, estimate, and error
df$SampleID = row.names(df)
df$Estimate=as.numeric(df$Estimate)
df$Error=as.numeric(df$Error)

# Merge the estimates with the sample data
RichPlot3 = merge(SamDat,df,by="SampleID")
head(RichPlot3)

## Option to save Breakaway's Richness Estimates
#RichPlot3$Model = as.character(RichPlot3$Model) # Model is a list, which prevents write.csv
#write.csv(RichPlot3,"Derived_data/CON18_Breakaway_Richness_Estimates.csv")

#### Load Breakaway Richness Estimates, if not continuing from above ####
#RichPlot3 = read.csv(file="Derived_data/CON18_Breakaway_Richness_Estimates.csv")


# Plot them a few ways
RichPlot3$TimesMixed = ordered(RichPlot3$TimesMixed, levels=c("Initial","1","2","4","8","16","32"))

p = ggplot(RichPlot3,aes(y=Estimate,x=TimesMixed,color=VortexControl))
p = p + geom_point()# + geom_errorbar(aes(ymin=Richness_estimate-Richness_stderr,ymax=Richness_estimate+Richness_stderr))
p

RichPlot = RichPlot3 %>%
    group_by(TimesMixed, VortexControl) %>%
    dplyr::summarize(Estimate_mean=mean(Estimate),Estimate_sd=sd(Estimate))

p = ggplot(RichPlot,aes(y=Estimate_mean,x=TimesMixed,color=VortexControl))
p = p + geom_point() + geom_errorbar(aes(ymin=Estimate_mean-Estimate_sd,ymax=Estimate_mean+Estimate_sd))
p = p + lims(y = c(0,3200))
p


#### Goal: summarize the estimates across different treatments (TimesMixed and VortexControl)

# First, make a single variable that has both of those elements
RichPlot3$Comp = paste(RichPlot3$TimesMixed,RichPlot3$VortexControl)
RichPlot3$Comp = as.factor(RichPlot3$Comp)
head(RichPlot3)

# Create an empty data frame that will hold the output data
RichPlotSumm = data.frame(Comp=levels(RichPlot3$Comp))
RichPlotSumm$Estimate = 0
RichPlotSumm$Error = 0
RichPlotSumm$p = 0
head(RichPlotSumm)

# Run a for loop that goes through each subset of data and makes the estimates using the betta function
for (i in levels(RichPlot3$Comp)){
  d = RichPlot3[RichPlot3$Comp==i,]
  Betta = betta(d$Estimate,d$Error)
  print(Betta$table)
  RichPlotSumm[RichPlotSumm$Comp==i,]$Estimate = Betta$table[,1]
  RichPlotSumm[RichPlotSumm$Comp==i,]$Error = Betta$table[,2]
  RichPlotSumm[RichPlotSumm$Comp==i,]$p = Betta$table[,3]
}

# Create a function to extract the rightmost values from our comparison variable (Comp)
# so we can get the TimesMixed and VortexControl info back
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n){
  substr(x, 1, nchar(x)-2)
}
# Pull out TimesMixed and VortexControl from the Comp variable so we can use them as plotting variables
RichPlotSumm$TimesMixed = substrLeft(as.character(RichPlotSumm$Comp),1)
RichPlotSumm$VortexControl = substrRight(as.character(RichPlotSumm$Comp),1)
RichPlotSumm
RichPlotSumm$GroupOrAgit = c("1", "16", "Agit", "2", "Agit", "32", "Agit", "4", "Agit", "8", "Agit", "Initial")

RichPlotSumm

#write.csv(RichPlotSumm,"Derived_data/CON18_Breakaway_Richness_Estimates_forplot.csv")

#### Load Breakaway Richness Estimates, if not continuing from above ####
#RichPlotSumm = read.csv(file="Derived_data/CON18_Breakaway_Richness_Estimates_forplot.csv")

RichPlotSumm$TimesMixed = ordered(RichPlotSumm$TimesMixed, levels=c("Initial","1","2","4","8","16","32"))
RichPlotSumm$GroupOrAgit = ordered(RichPlotSumm$GroupOrAgit, levels=c("Initial","1","2","4","8","16","32","Agit"))

# Plot final richness estimates with error bars
# ±1.96*SE represents 95% confidence intervals
p = ggplot(RichPlotSumm,aes(y=Estimate,x=TimesMixed,color=GroupOrAgit, shape=VortexControl))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylim(1000,3000) + ylab("Richness estimate")
p = p + xlab("Mixing frequency")
p = p + scale_x_discrete(labels=c("Initial"="Initial","1"="1\u00d7 (control)", "2"="2\u00d7", "4"="4\u00d7", "8"="8\u00d7", "16"="16\u00d7", "32"="32\u00d7"))
p = p + scale_color_manual(values = c("gray55", #initial
                                 "#440154FF", #1X
                                 "#404788FF",#2X
                                 "#2D708EFF",#4X 
                                 "#20A387FF",#8X 
                                 "#73D055FF",#16X 
                                 "#FDE725FF",#32X 
                                 "grey20"),#VortexControls
                               name="", guide=FALSE)
                               #labels=c("1X","2X","4X","8X","16X","32X","Vortex Control"))
p = p + scale_shape_manual(values=c(19,1), labels = c("Treatment (per X axis)", "Vortex control"),
                           guide = guide_legend(reverse = TRUE))
p = p + theme(legend.title = element_blank()) + theme(legend.position = c(0.25, 0.15), legend.box = "horizontal")
p # 475x350
