#### Get set up ####
library(corncob)
library(phyloseq)
library(dplyr)
library(ggplot2)

# Import phyloseq object
ps = readRDS("ps.CON18.remake")
ps

#### Differential abundance, SOIL MIXING TREATMENTS ####

# Create empty data frame to hold all results
df.mixed = data.frame()

# Subset for one comparison: a treatment group compared to baseline (1X)
# Loop to run through each times mixed
for (i in c("2","4","8","16","32")){
  # subset ps object to include 1x samples, plus all treatment samples
  ps.subset = prune_samples(sample_data(ps)$GroupOrAgit == "1" | sample_data(ps)$GroupOrAgit == i, ps)
  # Cut out any global zeros within this set
  ps.subset = prune_taxa(taxa_sums(ps.subset)>0,ps.subset)

  ## Drop taxa with super low abundance.
  # make a list of taxa that are relatively low abundance, then keep only the abundant taxa
  ps.subset.relabund = transform_sample_counts(ps.subset, function(x) x / sum(x))
  AbundTaxa = taxa_names(filter_taxa(ps.subset.relabund, function(x) mean(x) > 0.00001, TRUE))
  ps.subset = prune_taxa(AbundTaxa,ps.subset)
  
  # Create post-hoc filters from the relative abundance data
  # We ultimately won't be interested in enriched taxa that are rare even after enrichment
  ps.sub = prune_samples(sample_data(ps.subset.relabund)$GroupOrAgit==i,ps.subset.relabund)
  RareTrtTaxa = taxa_names(filter_taxa(ps.sub, function(x) mean(x) < 0.002, TRUE))
  
  # We give all parameters of interest (control and variable) to formula and phi.formula,
  # And then drop the parameter we want to test from the _null versions
  # (leaving 1 if there are no control variables, 
  # and the same parameters if we don't want to test for anything (as in phi.formula_null))
  # formula is the differential abundance
  # phi.formula is the differential variance
  # We may just need a very simple model for this dataset, testing for TimesMixed
  dT.ps.subset = differentialTest(formula = ~ TimesMixed, 
                               phi.formula = ~ TimesMixed,
                               formula_null = ~ 1,
                               phi.formula_null = ~ TimesMixed,
                               test = "Wald", boot = FALSE,
                               data = ps.subset,
                               fdr_cutoff = 0.05)
  
  # Making an empty dataframe to hold the full results
  df.ps.mixed  = data.frame()
  
  # Loop to pull out coefficients for each taxon
  for (j in 1:length(dT.ps.subset$significant_taxa)){
    # Get the significant model for that taxon
    sig_models = dT.ps.subset$significant_models[[j]]
    # Pull out the coefficients as above
    mu = data.frame(t(as.matrix(sig_models$coefficients[2,])))
    # Also grab the p_fdr estimate for that taxon's model
    p_fdr = dT.ps.subset$p_fdr[dT.ps.subset$significant_taxa][j]
    # Add that estimate onto our coefficient data frame
    mu$p_fdr = p_fdr
    # Create a column with the OTU ID
    mu$OTU= paste(row.names(data.frame(p_fdr)))
    # Add this row onto the df dataframe, which will collect the results
    # for all taxa as it iterates through this loop.
    df.ps.mixed = rbind(df.ps.mixed,mu)
  }
  
  # Clean up column names
  colnames(df.ps.mixed) = c("Estimate","SE","t","p","p_fdr","OTU")
  
  # Bring back in the taxonomy from the tax table
  SigOTUs = levels(as.factor(df.ps.mixed$OTU))
  pruned = prune_taxa(SigOTUs,ps.subset)
  taxtab = data.frame(tax_table(pruned))
  taxtab$OTU = c(taxa_names(pruned))
  joined.ps.subset = merge(df.ps.mixed,taxtab,by=c("OTU"))
  
  # Make column to designate if OTU is rare, for filtering later
  joined.ps.subset$RareTrtTaxa = ifelse(joined.ps.subset$OTU %in% RareTrtTaxa, "rare", "not rare")
  
  # Prep for merging with other times mixed
  joined.ps.subset$TimesMixed = paste(i)
  
  # Make final dataframe by joining together each differential test set
  df.mixed = rbind(df.mixed, joined.ps.subset)
}

dim(df.mixed)
head(df.mixed)

# Fix up some taxon naming issues
ignoreList = c("","uncultured","uncultured bacterium","uncultured soil bacterium","uncultured forest soil bacterium", "Ambiguous_taxa",
               "uncultured actinobacterium","uncultured planctomycete","uncultured Chloroflexi bacterium","uncultured Syntrophobacteraceae bacterium",
               "uncultured crenarchaeote","uncultured Gemmatimonadetes bacterium", "uncultured Acidobacteria bacterium",
               "uncultured Planctomyces sp.", "metagenome", "Subgroup 6", "Blastocatellia (Subgroup 4)", "uncultured Holophaga sp.",
               "uncultured Hyphomicrobiaceae bacterium", "uncultured proteobacterium", "1921-2", "RB41", "A21b",
               "AD3", "Subgroup_7", "Subgroup_2", "WD260")

df.mixed = df.mixed %>%
  mutate(Name = ifelse(Genus %in% ignoreList | is.na(Genus),
                       ifelse(Family %in% ignoreList | is.na(Family),
                              ifelse(Class %in% ignoreList | is.na(Class),
                                     ifelse(Phylum %in% ignoreList | is.na(Phylum),paste(OTU),
                                            paste(Phylum)),paste(Class)),paste(Family)),paste(Genus))) %>%
  mutate(Name = ifelse(Name == "uncultured thaumarchaeote","Thaumarchaeota",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus_Xiphinematobacter","Xiphinematobacter",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus Solibacter","Solibacter",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus_Udaeobacter","Udaeobacter",Name))

# Put treatments in order for graphing
df.mixed$TimesMixed = ordered(df.mixed$TimesMixed,levels=c("2", "4", "8", "16", "32"))

### Look at allll those OTUs
# How many total enriched (>0) (or depleted, <0) OTUs, total OR by treatment?
df.mixed$OTU = as.character(df.mixed$OTU)
total.OTUs=df.mixed%>%
  filter(Estimate<1)%>%
  count(TimesMixed) # for number of OTUs, by trt
  #summarize_each(funs(n_distinct)) # for total number of OTUs
total.OTUs

# How many big responders (mu > 1) that aren't overly rare?
total.OTUs=df.mixed%>%
  filter(Estimate>1)%>%
  filter(RareTrtTaxa=="not rare")%>%
  #count(TimesMixed) # for number of OTUs, by trt
  summarize_each(funs(n_distinct)) # for total number of OTUs
total.OTUs

# Tons of taxa; we will probably want to filter somewhat.
# Looking at biggest positive responders
df.mixed.BigResp.pos = df.mixed %>% filter(Estimate > 1)
dim(df.mixed.BigResp.pos)
# Remove the responders, who even after enrichment, are still 'rare'
df.mixed.BigResp.pos = df.mixed.BigResp.pos %>% filter(RareTrtTaxa=="not rare")
dim(df.mixed.BigResp.pos)
df.mixed.BigResp.pos

p = ggplot(df.mixed.BigResp.pos,aes(x=Estimate,color=Phylum,y= reorder(OTU, Estimate)))
p = p + theme_bw()
p = p + geom_point(size=2)
p = p + geom_errorbar(aes(xmax=Estimate+1.96*SE,xmin=Estimate-1.96*SE)) #plot negative estimate for y-value
p = p + xlab("Coefficient of differential abundance (vs. 1\u00d7 mixed)") + ylab("")
p = p + scale_y_discrete(breaks=df.mixed.BigResp.pos$OTU,labels=df.mixed.BigResp.pos$Name)
p = p + facet_grid(cols=vars(TimesMixed), scales="free_y", space="free_y")
p = p + theme(legend.position = "bottom")
p #750 x 650


# In order to make certain which treatment (1X or 32X) was used as the baseline, we can look for a certain taxa.
# We see Nocardioides as a responder, so it looks like positive estimate means enriched in 32X mixed.

#### Negative responders/Depleted taxa ####
# Retain only biggest negative responders
df.mixed.BigResp.neg = df.mixed %>% filter(Estimate < -1)
dim(df.mixed.BigResp.neg)

### Ultimately we are probably not interested in depleted taxa that were very rare to begin with...
# create a subset ps object for 1X communities, in relative abundance
ps.1X = prune_samples(sample_data(ps)$GroupOrAgit == "1", ps)
ps.1X = prune_taxa(taxa_sums(ps.1X)>0,ps.1X)
ps.1X = transform_sample_counts(ps.1X, function(x) x / sum(x))
# Create list of abundant 1X taxa.
Abundant1XTaxa = taxa_names(filter_taxa(ps.1X, function(x) mean(x) > 0.002, TRUE))

# Create sorting column based the OTU's abundance at 1X
df.mixed.BigResp.neg$Abund1XTaxa = ifelse(df.mixed.BigResp.neg$OTU %in% Abundant1XTaxa, "abundant1X", "rare1X")
head(df.mixed.BigResp.neg)

# keep only the negative responders that were abundant to begin with at 1X
df.mixed.BigResp.neg = df.mixed.BigResp.neg %>% filter(Abund1XTaxa=="abundant1X")
dim(df.mixed.BigResp.neg)
df.mixed.BigResp.neg

p = ggplot(df.mixed.BigResp.neg, aes(x=Estimate,color=Phylum,y= reorder(OTU, -Estimate)))
p = p + theme_bw()
p = p + geom_point(size=2)
p = p + geom_errorbar(aes(xmax=Estimate+1.96*SE,xmin=Estimate-1.96*SE)) #plot negative estimate for y-value
#p = p + theme(axis.text.x=element_text(angle=45,face="italic", size=6, hjust=1, vjust=1))
p = p + xlab("Coefficient of differential abundance (vs. 1\u00d7)") + ylab("")
p = p + scale_y_discrete(breaks=df.mixed.BigResp.neg$OTU,labels=df.mixed.BigResp.neg$Name)
p = p + facet_grid(cols=vars(TimesMixed), scales="free_y", space="free_y")
#p = p + geom_hline(yintercept=0)
p = p + theme(legend.position = "bottom")
p 

### So, we have the significant responders.
# Note - corncob gives us an estimate in a log-likelihood model that estimates relative abundance.
# "Estimate" is not "log2-fold change", but they are proportional.



#### Graph specific taxa of interest####

# create a ps for your treatment of interest, in relative abundance:
ps.subset = prune_samples(sample_data(ps)$GroupOrAgit == "1" | sample_data(ps)$GroupOrAgit == "32",ps)

ps.subset = prune_taxa(taxa_sums(ps.subset)>0,ps.subset)
ps.32.relabund = transform_sample_counts(ps.subset, function(x) x / sum(x))

# Prune the 32X relative abundance ps object down to just the Big Responders, and melt for graphing
ps.32.relabund.posresp = prune_taxa(df.mixed.BigResp.pos$OTU,ps.32.relabund)
ps.32.relabund.posresp = psmelt(ps.32.relabund.posresp)
head(ps.32.relabund.posresp)

#### Graph Positive responders  ####

# Create a list of OTUs for the 32X Big Responders
OTUNames=df.mixed.BigResp.pos%>%
  filter(TimesMixed=="32")%>%
  group_by(OTU,Phylum,Class,Order,Family,Genus)%>%
  summarize()
OTU32pos = c(paste(OTUNames$OTU))
OTU32pos

# take the almost full treatment dataset (1X and mixed communities)
# Cut out any global zeros within this set, transform to relative abundance
ps.fullish = prune_samples(sample_data(ps)$GroupOrAgit != "Agit" 
                           & sample_data(ps)$GroupOrAgit != "Initial",ps)
ps.fullish = prune_taxa(taxa_sums(ps.fullish)>0,ps.fullish)
ps.fullish.relabund = transform_sample_counts(ps.fullish, function(x) x / sum(x))
# Melt
ps.fullish.rel.melt = psmelt(ps.fullish.relabund)
ps.fullish.rel.melt


# Keep just the 32X pos responder OTUs in the almost full dataset
ps.fullish.rel.melt$OTU = as.factor(ps.fullish.rel.melt$OTU)
plot = ps.fullish.rel.melt %>%
  filter(OTU %in% OTU32pos)

# Figure of relative abundances (across all trts) for the 32X pos responders
# First, set up a Name column so we can include that in the facet labels...
# Because we have multiple OTUs that ultimately are only defined to, say, the genus level, these OTUs would
# be put in the same facet of the graph if there are multiple fromt eh same genus. 
# So we will create a faceting variable "both.names" which includes the OTU,
# and we will just insert a loooooong space so that we don't actually see the arbitrary OTU name in the facet label.
plot = plot %>%
  mutate(Name = ifelse(Genus %in% ignoreList | is.na(Genus),
                       ifelse(Family %in% ignoreList | is.na(Family),
                              ifelse(Class %in% ignoreList | is.na(Class),
                                     ifelse(Phylum %in% ignoreList | is.na(Phylum),paste(OTU),
                                            paste(Phylum)),paste(Class)),paste(Family)),paste(Genus)))

plot$both.names = paste0(as.character(plot$Name),"                          (", as.character(plot$OTU),")")

p = ggplot(plot)
p = p + geom_boxplot(aes(y=Abundance,x=TimesMixed, color=TimesMixed))
p = p + facet_wrap(~both.names, scales="free", nrow=3) +
  theme(strip.text = element_text(hjust=0)) +
  xlab("Mixing frequency") + ylab("Relative abundance")
p


# DF of all positive responding OTUs (all trts), with taxonomic info
OTUs=df.mixed.BigResp.pos%>%
  group_by(OTU,Phylum,Class, Order, Family,Genus,Species,Name, TimesMixed, Estimate, SE, t, p, p_fdr)%>%
  summarize()
OTUs
head(df.mixed.BigResp.pos)

#### Graph Negative responders ####
# There are a LOT of 32X responders, so we will graph the 16X negative responders.
# Create a list of OTUs for the 16X Big Responders, depleted
OTUNames=df.mixed.BigResp.neg%>%
  filter(TimesMixed=="32")%>%
  group_by(OTU,Phylum,Class,Order,Family,Genus)%>%
  summarize()
OTU32neg = c(paste(OTUNames$OTU))
OTU32neg

# take the almost full treatment dataset (1X and mixed communities)
# Cut out any global zeros within this set, transform to relative abundance
ps.fullish = prune_samples(sample_data(ps)$GroupOrAgit != "Agit" 
                           & sample_data(ps)$GroupOrAgit != "Initial",ps)
ps.fullish = prune_taxa(taxa_sums(ps.fullish)>0,ps.fullish)
ps.fullish.relabund = transform_sample_counts(ps.fullish, function(x) x / sum(x))
# Melt
ps.fullish.rel.melt = psmelt(ps.fullish.relabund)
ps.fullish.rel.melt


# Keep just the 32X pos responder OTUs in the almost full dataset
ps.fullish.rel.melt$OTU = as.factor(ps.fullish.rel.melt$OTU)
plot = ps.fullish.rel.melt %>%
  filter(OTU %in% OTU32neg)

# Figure of relative abundances (across all trts) for the 32X pos responders
# First, set up a Name column so we can include that in the facet labels
plot = plot %>%
  mutate(Name = ifelse(Genus %in% ignoreList | is.na(Genus),
                       ifelse(Family %in% ignoreList | is.na(Family),
                              ifelse(Class %in% ignoreList | is.na(Class),
                                     ifelse(Phylum %in% ignoreList | is.na(Phylum),paste(OTU),
                                            paste(Phylum)),paste(Class)),paste(Family)),paste(Genus)))%>%
  mutate(Name = ifelse(Name == "uncultured thaumarchaeote","Thaumarchaeota",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus_Xiphinematobacter","Xiphinematobacter",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus Solibacter","Solibacter",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus_Udaeobacter","Udaeobacter",Name))

plot$both.names = paste0(as.character(plot$Name),"                         ", as.character(plot$OTU),")")

p = ggplot(plot)
p = p + geom_boxplot(aes(y=Abundance,x=TimesMixed, color=TimesMixed))
p = p + facet_wrap(~both.names, scales="free", nrow=5) +
  theme(strip.text = element_text(hjust=0, size=7)) +
  xlab("Mixing frequency") + ylab("Relative abundance")
p



#### Differential abundance in the VORTEX CONTROLs ####
# This is largely a repeat of all the code above, tailored for the vortex control samples.
# [Agit is used to mean Vortex Controls]

# Create empty data frame to hold all results
df.vortex = data.frame()

# Subset for one comparison: a vortex treatment group compared to baseline (1X)
# Loop to run through each times mixed
for (i in c("2","4","8","16","32")){
  ps.subset = prune_samples(sample_data(ps)$GroupOrAgit == "1" | sample_data(ps)$GroupOrAgit == "Agit", ps) %>%
     prune_samples(sample_data(.)$GroupOrAgit == "1" | sample_data(.)$TimesMixed == i,.)
  # Cut out any global zeros within this set
  ps.subset = prune_taxa(taxa_sums(ps.subset)>0,ps.subset)
  # Drop taxa with super low abundance, and keep only the adundant taxa.
  ps.subset.relabund = transform_sample_counts(ps.subset, function(x) x / sum(x))
  AbundTaxa = taxa_names(filter_taxa(ps.subset.relabund, function(x) mean(x) > 0.00001, TRUE))
  ps.subset = prune_taxa(AbundTaxa,ps.subset)
  
  # Create post-hoc filters from the relative abundance data
  # We utimately won't be interested in enriched taxa that are rare even after enrichment
  ps.sub = prune_samples(sample_data(ps.subset.relabund)$GroupOrAgit == "Agit"
                           & sample_data(ps.subset.relabund)$TimesMixed == i, ps.subset.relabund)
  RareTrtTaxa = taxa_names(filter_taxa(ps.sub, function(x) mean(x) < 0.002, TRUE))
  
  # [see above version of the code for expanded explanations]
  dT.ps.subset = differentialTest(formula = ~ TimesMixed, 
                                  phi.formula = ~ TimesMixed,
                                  formula_null = ~ 1,
                                  phi.formula_null = ~ TimesMixed,
                                  test = "Wald", boot = FALSE,
                                  data = ps.subset,
                                  fdr_cutoff = 0.05)

  df.ps.vortex  = data.frame()
  # Loop to pull out coefficients for each taxon
  for (j in 1:length(dT.ps.subset$significant_taxa)){
    # Get the significant model for that taxon
    sig_models = dT.ps.subset$significant_models[[j]]
    mu = data.frame(t(as.matrix(sig_models$coefficients[2,])))
    p_fdr = dT.ps.subset$p_fdr[dT.ps.subset$significant_taxa][j]
    mu$p_fdr = p_fdr
    mu$OTU= paste(row.names(data.frame(p_fdr)))
    df.ps.vortex = rbind(df.ps.vortex,mu)
  }

  dim(df.ps.vortex)

  #Fix column names 
  colnames(df.ps.vortex) = c("Estimate","SE","t","p","p_fdr","OTU")
  
  # Let's bring back in our taxonomy from the tax table
  SigOTUs = levels(as.factor(df.ps.vortex$OTU))
  pruned = prune_taxa(SigOTUs,ps.subset)
  taxtab = data.frame(tax_table(pruned))
  taxtab$OTU = c(taxa_names(pruned))
  joined.ps.subset = merge(df.ps.vortex,taxtab,by=c("OTU"))
  
  # Make column to designate if OTU is rare, for filtering later
  joined.ps.subset$RareTrtTaxa = ifelse(joined.ps.subset$OTU %in% RareTrtTaxa, "rare", "not rare")
  
  # Prep for merging with other times mixed
  joined.ps.subset$TimesMixed = paste(i)
  
  # Make final dataframe by joining together each differential test set
  df.vortex = rbind(df.vortex, joined.ps.subset)
}

dim(df.vortex)
head(df.vortex)

# Fix up some taxon naming issues
ignoreList = c("","uncultured","uncultured bacterium","uncultured soil bacterium","uncultured forest soil bacterium", "Ambiguous_taxa",
               "uncultured actinobacterium","uncultured planctomycete","uncultured Chloroflexi bacterium","uncultured Syntrophobacteraceae bacterium",
               "uncultured crenarchaeote","uncultured Gemmatimonadetes bacterium", "uncultured Acidobacteria bacterium",
               "uncultured Planctomyces sp.", "metagenome", "Subgroup 6", "Blastocatellia (Subgroup 4)", "uncultured Holophaga sp.",
               "uncultured Hyphomicrobiaceae bacterium", "uncultured proteobacterium", "alphal_cluster", "1921-2", "JG30-KF-AS9")

df.vortex = df.vortex %>%
  mutate(Name = ifelse(Genus %in% ignoreList | is.na(Genus),
                       ifelse(Family %in% ignoreList | is.na(Family),
                              ifelse(Class %in% ignoreList | is.na(Class),
                                     ifelse(Phylum %in% ignoreList | is.na(Phylum),paste(OTU),
                                            paste(Phylum)),paste(Class)),paste(Family)),paste(Genus))) %>%
  mutate(Name = ifelse(Name == "uncultured_thaumarchaeote","Thaumarchaeota",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus_Xiphinematobacter","Xiphinematobacter",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus_Solibacter","Solibacter",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus_Ovatusbacter","Ovatusbacter",Name)) %>%
  mutate(Name = ifelse(Name == "Candidatus_Udaeobacter","Udaeobacter",Name))


# How many total enriched (>0) (or depleted, <0) OTUs, total OR by treatment?
total.OTUs=df.vortex%>%
  filter(Estimate>1)%>%
  count(TimesMixed) # for number of OTUs, by trt
  #summarize_each(funs(n_distinct)) # for total number of OTUs
total.OTUs

# How many big responders (>1) that aren't overly rare?
total.OTUs=df.vortex%>%
  filter(Estimate>1)%>%
  filter(RareTrtTaxa=="not rare")%>%
  #count(TimesMixed) # for number of OTUs, by trt
  summarize_each(funs(n_distinct)) # for total number of OTUs
total.OTUs

# Tons of taxa, will  want to filter somewhat.
# Looking at biggest positive responders
df.vortex.BigResp.pos = df.vortex %>% filter(Estimate > 1)
dim(df.vortex.BigResp.pos)

# Remove the responders, who even after enrichment, are still 'rare'
df.vortex.BigResp.pos = df.vortex.BigResp.pos %>% filter(RareTrtTaxa=="not rare")
dim(df.vortex.BigResp.pos)
df.vortex.BigResp.pos

p = ggplot(df.vortex.BigResp.pos,aes(x=Estimate,color=Phylum,y= reorder(OTU, Estimate)))
p = p + theme_bw()
p = p + geom_point(size=2)
p = p + geom_errorbar(aes(xmax=Estimate+1.96*SE,xmin=Estimate-1.96*SE)) #plot negative estimate for y-value
p = p + xlab("Coefficient of differential abundance (vs. 1\u00d7 mixed)") + ylab("")
p = p + scale_y_discrete(breaks=df.vortex.BigResp.pos$OTU,labels=df.vortex.BigResp.pos$Name)
p = p + facet_grid(cols=vars(TimesMixed), scales="free_y", space="free_y")
p = p + theme(legend.position = "bottom")
p

### Negative responders/Depleted taxa
# Looking at biggest negative responders
df.vortex.BigResp.neg = df.vortex %>% filter(Estimate < -1)
dim(df.vortex.BigResp.neg)

## Ultimately we are probably not interested in depleted taxa that were very rare to begin with...
## SO let's identify the originally rare taxa.

# create a subset ps object for 1X communities, in relative abundance
ps.1X = prune_samples(sample_data(ps)$GroupOrAgit == "1", ps)
ps.1X = prune_taxa(taxa_sums(ps.1X)>0,ps.1X)
ps.1X = transform_sample_counts(ps.1X, function(x) x / sum(x))
# Create list of abundant 1X taxa.
Abundant1XTaxa = taxa_names(filter_taxa(ps.1X, function(x) mean(x) > 0.002, TRUE))

# Create sorting column based the OTU's abundance at 1X
df.vortex.BigResp.neg$Abund1XTaxa = ifelse(df.vortex.BigResp.neg$OTU %in% Abundant1XTaxa, "abundant1X", "rare1X")
head(df.vortex.BigResp.neg)

# keep only the negative responders that were abundant to begin with at 1X
df.vortex.BigResp.neg = df.vortex.BigResp.neg %>% filter(Abund1XTaxa=="abundant1X")
dim(df.vortex.BigResp.neg)
df.vortex.BigResp.neg


p = ggplot(df.vortex.BigResp.neg, aes(x=Estimate,color=Phylum,y= reorder(OTU, Estimate)))
p = p + theme_bw()
p = p + geom_point(size=2)
p = p + geom_errorbar(aes(xmax=Estimate+1.96*SE,xmin=Estimate-1.96*SE)) #plot negative estimate for y-value
p = p + xlab("Coefficient of differential abundance (vs. 1\u00d7 mixed)") + ylab("")
p = p + scale_y_discrete(breaks=df.vortex.BigResp.neg$OTU,labels=df.vortex.BigResp.neg$Name)
p = p + facet_grid(cols=vars(TimesMixed), scales="free_y", space="free_y")
p = p + theme(legend.position = "bottom")
p 

### So, we have a list of the significant responders
# Note - corncob gives us an estimate in a log-likelihood model that estimates relative abundance.
# "Estimate" is not "log2-fold change", but they are proportional.

### Compare relative abundances of specific taxa of interest

# create a ps for your vortex treatment of interest, in relative abundance:
ps.subset.agit = prune_samples(sample_data(ps)$GroupOrAgit == "1" | sample_data(ps)$GroupOrAgit == "Agit", ps) %>%
  prune_samples(sample_data(.)$GroupOrAgit == "1" | sample_data(.)$TimesMixed == "32",.)

ps.subset.agit = prune_taxa(taxa_sums(ps.subset.agit)>0,ps.subset.agit)
ps.32.relabund = transform_sample_counts(ps.subset.agit, function(x) x / sum(x))

# Prune the 32X relative abundance ps object down to just the Big Responders, and melt for graphing
ps.32.relabund.posresp.vort = prune_taxa(df.vortex.BigResp.pos$OTU,ps.32.relabund)
ps.32.relabund.posresp.vort = psmelt(ps.32.relabund.posresp.vort)
head(ps.32.relabund.posresp.vort)

# Graph comparing RelAbund of 32X Big Responders, 1X and 32X
p = ggplot(ps.32.relabund.posresp.vort)
p = p + geom_boxplot(aes(y=Abundance,x=TimesMixed, color=TimesMixed))
p = p + facet_wrap(~OTU,scales="free")
p

### graph relative abundance of the 32X pos responders across all treatments
# take the almost full treatment dataset (1X and mixed communities)
ps.fullish = prune_samples(sample_data(ps)$GroupOrAgit != "Agit" 
                           & sample_data(ps)$GroupOrAgit != "Initial",ps)
ps.fullish
# Cut out any global zeros within this set
ps.fullish = prune_taxa(taxa_sums(ps.fullish)>0,ps.fullish)
ps.fullish
# Transform to relative abundance
ps.fullish.relabund = transform_sample_counts(ps.fullish, function(x) x / sum(x))
# Melt
ps.fullish.rel.melt = psmelt(ps.fullish.relabund)
ps.fullish.rel.melt

# Create a list of OTUs for the 32X Big Responders
OTUNames.vortex=ps.32.relabund.posresp.vort%>%
  group_by(OTU,Phylum,Class,Order,Family,Genus)%>%
  summarize()
OTUNames.vortex
OTU32pos.vortex = c(paste(OTUNames.vortex$OTU))

# Keep just the 32X pos responder OTUs in the almost full dataset
ps.fullish.rel.melt$OTU = as.factor(ps.fullish.rel.melt$OTU)
plot.vortex = ps.fullish.rel.melt %>%
  filter(OTU %in% OTU32pos.vortex)

plot.vortex = plot.vortex %>%
  mutate(Name = ifelse(Genus %in% ignoreList | is.na(Genus),
                       ifelse(Family %in% ignoreList | is.na(Family),
                              ifelse(Class %in% ignoreList | is.na(Class),
                                     ifelse(Phylum %in% ignoreList | is.na(Phylum),paste(OTU),
                                            paste(Phylum)),paste(Class)),paste(Family)),paste(Genus)))#%>%

plot$both.names = paste0(as.character(plot$Name),"                          (", as.character(plot$OTU),")")

p = ggplot(plot.vortex)
p = p + geom_boxplot(aes(y=Abundance,x=TimesMixed, color=TimesMixed))
p = p + facet_wrap(~both.names,scales="free", nrow=2) +
  theme(strip.text = element_text(hjust=0))+
  xlab("Mixing frequency") + ylab("Relative abundance")
p


OTUs.vortex=plot.vortex%>%
  group_by(OTU,Phylum,Family,Genus,Name)%>%
  summarize()
OTUs.vortex
