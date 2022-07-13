# Load packages
library(phyloseq)
library(vegan)
library(ggplot2)
library(magrittr)
library(DescTools) # for Dunnett's test

# Bring in ps object
ps.TRT = readRDS("ps.CON18.remake")
# Normalize to relative abundances.
ps.TRT.norm = transform_sample_counts(ps.TRT, function(x) x / sum(x) )
# rename for ease of use.
ps = ps.TRT.norm

#### PCoA ####

# Create PCoA ordination
ps.ordination.PCoA = ordinate(ps, method="PCoA", distance="bray")

# Order treatments/mixing sets (combined variable "Group2") for graphing
sample_data(ps)$Group2 = ordered(sample_data(ps)$Group2, 
                                 levels=c("Initial","1X","U","V","W","X","Q","R","S","T","M",
                                          "N","O","P","I","J","K","L","E","F","G","H",
                                          "ZAgit2","ZAgit4","ZAgit8","ZAgit16","ZAgit32"))
# plot ordination
p = plot_ordination(ps, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "Group2", shape = "VortexControl", label = NULL, justDF = FALSE)
p = p + theme_bw() + geom_point(size =1.5)
p = p + scale_color_manual(breaks = c("Initial", "1X","U","Q","M","I","E"),
                           values = c("gray65", #initial
                                      "#440154FF", #1X
                                      "#404788FF","#404788FF","#404788FF","#404788FF",#2X
                                      "#2D708EFF","#2D708EFF","#2D708EFF","#2D708EFF",#4X 
                                      "#20A387FF","#20A387FF","#20A387FF","#20A387FF",#8X 
                                      "#73D055FF","#73D055FF","#73D055FF","#73D055FF",#16X 
                                      "#FDE725FF","#FDE725FF","#FDE725FF","#FDE725FF",#32X 
                                      "#404788FF","#2D708EFF","#20A387FF","#73D055FF","#FDE725FF"),#VortexControls
                           name="", 
                           labels=c("Initial","1\u00d7 (control)","2\u00d7","4\u00d7","8\u00d7","16\u00d7","32\u00d7"))
p = p + scale_shape_manual(breaks = c("Y", "N"), name = "",
                           values = c(1, 19), labels = c("Vortex control", ""))
p = p + theme(legend.text = element_text(size = 10))
p$layers = p$layers[-1] #to remove the larger point coded in original plot
p

## PERMANOVA on BC-Dissimilaity (PCoA)
# We want to know whether the distances between samples correspond to the treatments

# Create a variable that is your distance matrix; You'll need this later.
# (exclude the initial community and vortex controls)
ps.trt = subset_samples(ps, TimesMixed != "Initial" & VortexControl=="N")
DistVar = phyloseq::distance(ps.trt, method = "bray")

# Extract the sample_data from the phyloseq object and turn it into a dataframe; You'll need this later.
psData = data.frame(sample_data(ps.trt))

# Run PERMANOVA by TimesMixed
adonis(DistVar ~ TimesMixed, data = psData, method = "bray")

# Create a function for pairwise comparisons, to test significant effects between each treatment
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='BH',reduce=NULL,perm=999)
{
  co <- combn(unique(as.character(factors)),2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  for(elem in 1:ncol(co)){
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    }
    
    else  (
      if (sim.function == 'daisy'){
        x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      } 
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    ad <- adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])],
                 permutations = perm);
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad$aov.tab[1,1])
    SumsOfSqs <- c(SumsOfSqs, ad$aov.tab[1,2])
    F.Model <- c(F.Model,ad$aov.tab[1,4]);
    R2 <- c(R2,ad$aov.tab[1,5]);
    p.value <- c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
  
  if(!is.null(reduce)){
    pairw.res <- subset (pairw.res, grepl(reduce,pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
    
    sig = c(rep('',length(pairw.res$p.adjusted)))
    sig[pairw.res$p.adjusted <= 0.1] <-'.'
    sig[pairw.res$p.adjusted <= 0.05] <-'*'
    sig[pairw.res$p.adjusted <= 0.01] <-'**'
    sig[pairw.res$p.adjusted <= 0.001] <-'***'
    pairw.res <- data.frame(pairw.res[,1:7],sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
} 

## Run pairwise adonis function on distance variable, if there was a significant overall effect
pairwise.adonis(DistVar,psData$TimesMixed)


#### Bray-Curtis dissimilarities: mixed communities compared to 1X communities ####
# Create function to expand.grid without duplicate pairs (such as a,b and b,a) or equals (i.e. a,a)
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  y <- unique(y)
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

# Expand grid for comparisons of interest: everything vs 1X communities
X1IDs = sample_names(subset_samples(ps,TimesMixed=="1"))
X1Comp = sample_names(ps)
x = X1IDs
y = X1Comp
BC = expand.grid.unique(x,y)
head(BC)
colnames(BC) = c("X1Comp","SampleID")
BC = as.data.frame(BC)
BC$X1Comp = as.factor(BC$X1Comp)
BC$SampleID = as.factor(BC$SampleID)

# Find BC Dissimilarity (this takes several minutes++)
for (i in 1:nrow(BC)){
  ps.BC = subset_samples(ps,SampleID == BC$X1Comp[i] | SampleID ==BC$SampleID[i])
  ps.BC = prune_taxa(taxa_sums(ps.BC)>0,ps.BC)
  BC$Dissimilarity[i] = vegdist(otu_table(ps.BC),method="bray")[1]
}
head(BC)

#### Graph of Bray-Curtis Dissimilarities compared to 1X communities ####
# Bring in sample data, and order the plotting factors
BC.plot = merge(data.frame(BC),data.frame(sample_data(ps)),by="SampleID")
BC.plot$TimesMixed = ordered(BC.plot$TimesMixed, levels=c("Initial","1","2","4","8","16","32"))

p.0 = ggplot(subset(BC.plot, TimesMixed != "Initial"),aes(x=TimesMixed,y=Dissimilarity,color=GroupOrAgit))
p.0 = p.0 + geom_jitter(alpha = 0.025)
p.0 = p.0 + geom_boxplot() + expand_limits(y=c(0.05,0.85)) +
  labs(x="Mixing frequency", y="Dissimilarity, compared to 1\u00d7", title=NULL)
p.0 = p.0 + scale_x_discrete(labels=c("1"="1\u00d7 (control)", "2"="2\u00d7", "4"="4\u00d7", "8"="8\u00d7", "16"="16\u00d7", "32"="32\u00d7"))
p.0 = p.0 + scale_color_manual(values = c(
  "#440154FF", #1X
    "#404788FF",#2X
    "#2D708EFF",#4X 
    "#20A387FF",#8X 
    "#73D055FF",#16X 
    "#FDE725FF",#32X 
    "grey20"),#VortexControls
  name="", 
  labels=c("1\u00d7 (control)","2\u00d7","4\u00d7","8\u00d7","16\u00d7","32\u00d7","Vortex control"))
p.0 = p.0 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p.0 = p.0 + theme(legend.position = "none")
#p.0 = p.0 + theme(legend.position=c(.75,.78))
p.0 #475x275

# Statistics for comparison vs 1X are above (PERMANOVA and pairwise adonis)

#### Bray-Curtis dissimilarities, within mixing sets ####
# You'll be using the same function, expand.grid.unique, created above.

# Expand grid for comparisons of interest: within the same "Group2" factor
# i.e., within the same mixing set, or within 1X, or within the same vortex control trt
sampnames = sample_names(ps)
x = sampnames
y = sampnames
BC.InSets = expand.grid.unique(sampnames, sampnames)

# Clean things up
colnames(BC.InSets) = c("sample1","sample2")
BC.InSets = as.data.frame(BC.InSets)
BC.InSets$sample1 = as.factor(BC.InSets$sample1)
BC.InSets$sample2 = as.factor(BC.InSets$sample2)
dim(BC.InSets)

# Pull in mixing set numbers for sample1, then sample2
SamData = data.frame(sample_data(ps)[,c("SampleID","Group2")])
colnames(SamData)=c("sample1","sample1.set")
BC.InSets = merge(BC.InSets,SamData,by="sample1")

SamData = data.frame(sample_data(ps)[,c("SampleID","Group2")])
colnames(SamData)=c("sample2","sample2.set")
BC.InSets = merge(BC.InSets,SamData,by="sample2")
head(BC.InSets)

str(BC.InSets) # sample numbers need to be characters, and sets need to be factors
BC.InSets$sample2 = as.character(BC.InSets$sample2)
BC.InSets$sample1 = as.character(BC.InSets$sample1)

# Make sorting columns to determine is samples are in the same mixing set, and to rule out same sample comparisons
BC.InSets$sameset = ifelse(BC.InSets$sample1.set==BC.InSets$sample2.set, "same", "no")
BC.InSets$samesample = ifelse(BC.InSets$sample1==BC.InSets$sample2, "same", "no")

# Create a new df with just your comparisons of interest
BC.InSets2 = BC.InSets %>%
  dplyr::filter(sameset == "same" & samesample == "no")
dim(BC.InSets2)

# Calculate Bray-Curtis dissimilarities for your within-mixing set comparisons (this may take >10 mins)
for (i in 1:nrow(BC.InSets2)){
  ps.BC = subset_samples(ps,SampleID == BC.InSets2$sample1[i] | SampleID ==BC.InSets2$sample2[i])
  ps.BC = prune_taxa(taxa_sums(ps.BC)>0,ps.BC)
  BC.InSets2$Dissimilarity[i] = vegdist(otu_table(ps.BC),method="bray")[1]
}

head(BC.InSets2)
#### Graph of Bray-Curtis-Dissimilarities, within Mixing Sets (and within 1X, vortex controls by trt)
# Add sample data, and order factors for graphing
BC.InSets2$SampleID = BC.InSets2$sample2
BC.InSets.plot = merge(data.frame(BC.InSets2),data.frame(sample_data(ps)),by="SampleID")
BC.InSets.plot$TimesMixed = ordered(BC.InSets.plot$TimesMixed, levels=c("Initial","1","2","4","8","16","32"))

# Create plot
p.3 = ggplot(subset(BC.InSets.plot, TimesMixed != "Initial"),aes(x=TimesMixed,y=Dissimilarity,color=GroupOrAgit))
p.3 = p.3 + geom_jitter(alpha = 0.05)
p.3 = p.3 + geom_boxplot() + expand_limits(y=c(0.05,0.85)) +
  labs(x="Mixing frequency", y="Dissimilarity, within pooled mixing set", title=NULL)
p.3 = p.3 + scale_x_discrete(labels=c("1"="1\u00d7 (control)", "2"="2\u00d7", "4"="4\u00d7", "8"="8\u00d7", "16"="16\u00d7", "32"="32\u00d7"))
p.3 = p.3 + scale_color_manual(values = c(
                                 "#440154FF", #1X
                                 "#404788FF",#2X
                                 "#2D708EFF",#4X 
                                 "#20A387FF",#8X 
                                 "#73D055FF",#16X 
                                 "#FDE725FF",#32X 
                                 "grey20"),#VortexControls
                               name="", 
                               labels=c("1\u00d7 (control)","2\u00d7","4\u00d7","8\u00d7","16\u00d7","32\u00d7","Vortex control"))
p.3 = p.3 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(strip.background = element_rect(fill="white"))
p.3 = p.3 + theme(legend.position=c(.75,.75)) + theme(legend.text=element_text(size=8)) +
  theme(legend.margin=margin(t=-0.5, r=0, b=0, l=0, unit="cm"))
p.3 #475x275

# Statistical comparison of disperson within each pooled mixing set
# using betadisper, which is more robust than ANOVA of BC-Dis values
# See Anderson (2006) for a discussion on tests of multivariate dispersion.
# Run this separately for each treatment
beta = betadisper(DistVar, psData$Group)
anova(beta)
test = permutest(beta, pairwise=TRUE, permutations=999)
test
Tukey=TukeyHSD(beta, conf.level = 0.95)
Tukey

