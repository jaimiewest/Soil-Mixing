#### Getting Set Up ####
library(BiocManager)
library(phyloseq)
library(vegan)
library(magrittr)
library(dplyr)
library(ggplot2)

# Bring in ps object
ps.TRT = readRDS("ps.CON18.remake")

# Normalize to relative abundances. 
ps.TRT.NP = transform_sample_counts(ps.TRT, function(x) x / sum(x) )

# Create a melted df for use below
df.melt = psmelt(ps.TRT.NP) 

#### Determine the most abundant phyla ####
# Get abundant phyla names by:Group by each sample, and phylum within that, and sum up (previously normalized) counts.
# Then, across all samples together (group_by(Phylum)), get mean phyla abundances across all samples. And sort/arrange.
Phy = df.melt%>%
  #filter(TimesMixed == "32" & VortexControl == "N")%>%
  filter(TimesMixed == "Initial")%>%
  group_by(Sample, Phylum)%>%
  dplyr::summarize(PhySum = sum(Abundance))%>%
  group_by(Phylum)%>%
  dplyr::summarize(PhyMean = mean(PhySum))%>%
  arrange(-PhyMean)
Phy

# Select phyla to keep; try the top 10 for each trt.
PhytoKeepInit = c(paste(Phy$Phylum[c(1:10)]))

PhytoKeep = c(PhytoKeepInit, PhytoKeep1, PhytoKeep2, PhytoKeep4, PhytoKeep8, PhytoKeep16, PhytoKeep32)
PhytoKeep = unique(PhytoKeep)
PhytoKeep
# Aiming for 12 so let's remove WPS-2
PhytoKeep = PhytoKeep[c(1:11,13)]
PhytoKeep

#### Box plots, facet by PHYLA--SI FIGURE 4 ####
df.plot = df.melt %>%
  group_by(Sample, Phylum, TimesMixed, VortexControl) %>% 
  dplyr::summarize(PhySum=sum(Abundance)) %>%
  filter(Phylum %in% PhytoKeep) 

p2 = ggplot(df.plot,aes(x=TimesMixed,y=PhySum,color=VortexControl))
p2 = p2 + geom_boxplot()
p2 = p2 + facet_wrap(~Phylum, scales="free")
p2 = p2 + theme(strip.text = element_text(size=10,face="italic")) +
  scale_x_discrete(labels=c("Initial"="Initial", "1"="1\u00d7",
                            "2"="2\u00d7", "4"="4\u00d7", "8"="8\u00d7", "16"="16\u00d7", "32"="32\u00d7")) +
  labs(x = "Mixing frequency", y = "Relative abundance", fill = "") +
  scale_color_manual(values = c("darkslategray",
                                  "darkseagreen3"),
                       name="", 
                       labels=c("Mixed soil treatment","Vortex control")) +
  theme_bw()
p2


##### Make a list of the most abundant OTUs in 32X ####
OTU32 = df.melt%>%
  filter(TimesMixed == "32" & VortexControl == "N")%>%
  group_by(Sample, OTU)%>%
  dplyr::summarize(OTUSum = sum(Abundance))%>%
  arrange(-OTUSum)
OTU32



# this snipet shows how the top 20 OTUs at 1X make up 25% of mean relative abundance,
# whereas the top 20 OTUs at 32X make up almost 70% of mean relative abundance at 32x. 
# (see rank-abundance Jupyter notebook for full rank abundance curves)

TopTaxa.1x = df.melt%>%
  filter(TimesMixed == "1" & VortexControl == "N")%>%
  filter(Abundance > 0)%>%
  group_by(OTU)%>%
  dplyr::summarize(MeanRelAb = mean(Abundance))%>%
  arrange(-MeanRelAb)
sum(TopTaxa.1x$MeanRelAb[1:20])

TopTaxa.32x = df.melt%>%
  filter(TimesMixed == "32" & VortexControl == "N")%>%
  filter(Abundance > 0)%>%
  group_by(OTU)%>%
  dplyr::summarize(MeanRelAb = mean(Abundance))%>%
  arrange(-MeanRelAb)
sum(TopTaxa.32x$MeanRelAb[1:20])
