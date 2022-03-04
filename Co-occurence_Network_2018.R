# Creating networks to identify co-occuring bacteria (and co-exclusions)
# Modified from code used in Whitman et al. (2019)
# See Connor et al. (2017) for rationale behind approach to network construction

# Load packages
library(phyloseq)
library(ggplot2)
library(reshape)
library(plyr)
library(dplyr)
library(vegan)
library(ade4)
library(wesanderson)
library(igraph)
library(tidyr)
library(threejs)
library(htmlwidgets)
library(RColorBrewer)

#### Bringing in data ####
# Bring in ps object
ps.TRT = readRDS("ps.CON18.remake")
# transform to relative abundance
ps.TRT.norm = transform_sample_counts(ps.TRT, function(x) x/sum(x))

# We want to remove the least abundant taxa, to reduce computational load and to not bother with low-abundance taxa
# We could interpret this as being, across all samples, at least X abundant
# That way, something low abundance but widely present, or very high abundance, will get included
# The cutoff is, of course, somewhat arbitrary - present across all samples at 0.005 seems quite low, though.
cutoff = 0.005

ps.mini = prune_taxa((taxa_sums(ps.TRT.norm) > cutoff), ps.TRT.norm)

# May want to trim ps object (e.g., if making network from 4x samples only)
ps.trim = subset_samples(ps.mini,TimesMixed == "4")
ps.trim = subset_samples(ps.trim, VortexControl == "N")
ps.trim = prune_taxa((taxa_sums(ps.trim) > 0), ps.trim)

# Set phyloseq object
ps=ps.trim

#### Create matrix with no noise added ####

# Creating a cutoff function to collect the Spearman rho cutoff 
# and the fraction of OTUs that are included in largest cluster
cutoff_function = function(cutoff,ps){
  adjacency_matrix = as.matrix(ps)
  adjacency_matrix[abs(adjacency_matrix)<cutoff] = 0
  adjacency_matrix[abs(adjacency_matrix)>cutoff] = 1
  am = graph.adjacency(adjacency_matrix)
  c = clusters(am)
  MaxSize = c$csize[1]
  return(c(cutoff,MaxSize/dim(adjacency_matrix)[1]))
}

# Test it against this range of values
# I first started out with wide range from 0 to 1, then narrowed it down to pinpoint optimal rho cutoff value
inputs = seq(0.4,0.8, by = 0.01)

# Setting up the data frame
Rep=c()
cutoff=c()
PercentIncl=c()
df.no=data.frame(Rep,cutoff,PercentIncl)

network_function = function(ps){
  for (i in 1:16){
    # Record iteration
    Rep = i
    # Calculate spearman correlations
    ps.dist = cor(otu_table(ps), use="everything", method="spearman") 
    # Apply the function to the range of values and turn it into a dataframe
    df = t(sapply(inputs, cutoff_function, ps=ps.dist))
    df = data.frame(df)
    df = data.frame(i,df)
    colnames(df) = colnames(df.no)
    df.no=rbind(df.no,df)
  }
  colnames(df.no) = c("Rep","Cutoff","PctIncl")
  return(df.no)
}

df.n = network_function(ps)

# Plot percent of OTUs included in largest cluster with increasing Spearman cutoffs
p = ggplot(df.n)
p = p + geom_point(aes(x=Cutoff,y=PctIncl, color=Rep))
p

# Dropoff around 0.61 for "Initial"
# Dropoff around 0.75 for 32X

#### Matrix with just noise added ####

## To calculate minimum distance between taxa across the whole dataset
# Basically, it will be 1/[the most total sequences across all samples]
mindist = min(1/max(sample_sums(ps)))
# Across all datasets, the minimum distance between taxa within a sample
delta = mindist*10^-2
# Just to be generous

# Setting up the data frame
Rep=c()
cutoff=c()
PercentIncl=c()
df.list=data.frame(Rep,cutoff,PercentIncl)

noise_function = function(ps){
  # Running it multiple times
  for (i in 1:16){
    # Record iteration
    Rep = i
    # Making matrix of random values that are less than smallest difference between two samples
    b = delta/1000
    E = replicate(dim(otu_table(ps))[1], rnorm(dim(otu_table(ps))[2]))
    E = 2*b*E
    E = -b + E
    # Add the noise to the matrix
    otu_table(ps) = otu_table(as.matrix(otu_table(ps))+abs(t(E)),taxa_are_rows=FALSE)
    # Calculate spearman correlations
    ps.dist = cor(otu_table(ps), use="everything", method="spearman")
    # Apply the function to the range of values and turn it into a dataframe
    df = t(sapply(inputs,cutoff_function,ps=ps.dist))
    df = data.frame(df)
    df = data.frame(i,df)
    colnames(df) = colnames(df.list)
    df.list=rbind(df.list,df)
  }
  colnames(df.list) = c("Rep","Cutoff","PctIncl")
  return(df.list)
}

df.e = noise_function(ps)

p = ggplot(df.e,aes(x=Cutoff,y=PctIncl, color=Rep))
p = p + geom_point()
p

# Now the dropoff starts a bit earlier; including random error to break ties 
# reduces the number of taxa included in the network.

#### Matrix with noise and permuations added ####

# Setting up the data frame
Rep=c()
Cutoff=c()
PctIncl=c()
df.list.permute=data.frame(Rep,Cutoff,PctIncl)

permute_noise_function = function(ps){
  # Running it multiple times
  for (i in 1:16){
    # Record iteration
    Rep = i
    # Making matrix of random values as above
    b = delta/1000
    E = replicate(dim(otu_table(ps))[1], rnorm(dim(otu_table(ps))[2]))
    E = 2*b*E
    E = -b + E
    
    # Add the noise to the matrix
    otu_table(ps) = otu_table(as.matrix(otu_table(ps))+abs(t(E)),taxa_are_rows=FALSE)
    # Get the OTU table as matrix
    M = matrix(otu_table(ps))
    # number of elements in matrix
    n = dim(M)[1]*dim(M)[2]
    # Sample the matrix with replacement; make new random matrix
    M.new = matrix(base::sample(M, n, replace=TRUE),nrow=dim(otu_table(ps))[1])
    # Calculate spearman correlations
    ps.dist = cor(M.new, use="everything", method="spearman") 
    # Apply the function to the range of values and turn it into a dataframe
    df = t(sapply(inputs, cutoff_function, ps=ps.dist))
    df = data.frame(df)
    df = data.frame(i,df)
    colnames(df) = colnames(df.list.permute)
    df.list.permute=rbind(df.list.permute,df)
  }
  colnames(df.list.permute) = c("Rep","Cutoff","PctIncl")
  return(df.list.permute)
}

df.p = permute_noise_function(ps)

p = ggplot(df.p,aes(x=Cutoff,y=PctIncl, color=Rep))
p = p + geom_point()
p

# As we can see, permuting the matrix randomly still makes networks.
# Next, we must choose a rho cutoff that is above the random threshold.

# Join the three types of runs together
df.n$Set = "Nothing"
df.p$Set = "Error+Permute"
df.e$Set = "Error"
df.full = rbind(df.n,df.e,df.p)

#df.full

# Try different values for cutoff. For 2018 dataset, I used:
# 0.61 for "Initial"
# 0.6 for 32X
# 0.61 for 16X
# 0.61 for 8X
# 0.61 for 4X
# 0.62 for 2X
# 0.62 for 1X
Cutoff=0.61

# Choose a cutoff value that keeps the max below 0.01 (1%) in the error+permuted dataset
max(df.full[df.full$Cutoff==Cutoff  & df.full$Set=="Error+Permute",]$PctIncl)

p = ggplot(df.full,aes(x=Cutoff,y=PctIncl, color=Set)) + geom_point() #+ ylim(c(0,0.1))
p

# The above code basically does what we want - calculates how large the largest cluster is 
# for each of the rho cutoff values, with no alterations, then adding random error, and finally adding 
# random error plus rearranging the matrix randomly to result in no expected correlations. 
# Thus, we can see where the rho cutoff values fall for all options. We should choose a cutoff 
# that is higher than where the permuted one reaches ~1% of the OTUs being included in the 
# largest cluster. We can then look at that cutoff in the actual (with added random error) data.

# To be really confident in this, any chosen cutoff should be examined to make sure 
# small variations in its value do not affect ecological conclusions

# Now we should be ready to use that cutoff to calculate the various parameters of the system.
# Now, to test our confidence in the overall network developed, we will caculate a series of 
# parameters about the network, and determine their null distributions.

# First, we will use the derived network, based on the cutoff that we chose above.

#### Matrix with just noise added ####

#### Create Network ####
# Set rho cutoff based on above calculations, if not continuing from above
Cutoff=0.61
rho = Cutoff

## set delta, if not continuing from above
mindist = min(1/max(sample_sums(ps))) # To calculate minimum distance between taxa across the whole dataset
delta = mindist*10^-2 # just to be generous

matrix_parameters = function(ps,rho){
  df.reps=data.frame(m=c(),n=c(),k=c(),apl=c(),c=c())
  # Ultimately, run 1000 times (maybe test it with just 100)
  for (i in 1:1000){
    # Record iteration
    Rep = i
    # Making matrix of random values
    b = delta/1000
    E = replicate(dim(otu_table(ps))[1], rnorm(dim(otu_table(ps))[2]))
    E = 2*b*E
    E = -b + E
    # Add the noise to the matrix
    otu_table(ps) = otu_table(as.matrix(otu_table(ps))+abs(t(E)),taxa_are_rows=FALSE)
    # Calculate spearman correlations
    ps.dist = cor(otu_table(ps), use="everything", method="spearman") 
    # Turns Spearman into matrix
    adjacency_matrix = as.matrix(ps.dist)
    # Set correlations below rho cutoff to 0, and above to 1
    adjacency_matrix[abs(adjacency_matrix)<rho] = 0
    adjacency_matrix[abs(adjacency_matrix)>rho] = 1
    # Remove all taxa that have no correlations above the cutoff.
    adjacency_matrix = adjacency_matrix[rowSums(adjacency_matrix)>0,]
    # Create adjacency matrix from correlation matrix
    adjacency_matrix = graph.adjacency(adjacency_matrix, diag=FALSE, mode=c("undirected"), weighted=TRUE)
    # Gets edge list from adjacency matrix
    edge_list = get.edgelist(adjacency_matrix)
    # Create the network
    N = graph_from_data_frame(edge_list,directed=FALSE)
    m = ecount(N) # Number of edges
    n = vcount(N) # Number of nodes
    k = 2*m/n # Average degree
    apl = mean_distance(N,directed=FALSE) # Average path length
    c = transitivity(N) # Clustering coefficient of whole graph
    
    # Apply the function to the range of values and turn it into a dataframe
    df = data.frame(m,n,k,apl,c)
    colnames(df) = colnames(df.reps)
    df.reps=rbind(df.reps,df)
  }
  colnames(df.reps)=c("m","n","k","apl","c")
  return(df.reps)
}

df.reps = matrix_parameters(ps,rho)

k.ave = mean(df.reps$k)
n.ave = mean(df.reps$n)
edges = k.ave * n.ave / 2
p = edges * 2 / (n.ave*n.ave) # Probability that any pair of verticies is connected 

# Thus we can generate a random matrix that preserves the degree of the network
ER_function = function(k.ave,n.ave,edges,p){
  df.ER=data.frame(Rep=c(),m=c(),n=c(),k=c(),apl=c(),c=c())
  for (i in 1:1000){
    Rep=i
    # Makes a random uniform matrix with the average number of nodes
    ER.thresh = matrix(runif(round(n.ave)*round(n.ave)),ncol=round(n.ave))
    # Set anything in the upper triangle that is below the threshold (1-p) to zero and anything above to 1
    ER.thresh[upper.tri(ER.thresh,diag=FALSE) & ER.thresh < (1-p)] = 0
    ER.thresh[upper.tri(ER.thresh,diag=FALSE) & ER.thresh > (1-p)] = 1
    # Reflect the upper triangle to the lower triangle
    ER.thresh[lower.tri(ER.thresh,diag=FALSE)] = 0
    # Set the diagonal to 0 (no self-connections)
    diag(ER.thresh) = 0
    # This is the ER random thresholded matrix
    
    # Create adjacency matrix from correlation matrix
    adjacency_matrix = graph.adjacency(ER.thresh, diag=FALSE, mode=c("undirected"), weighted=TRUE)
    # Gets edge list from adjacency matrix
    edge_list = get.edgelist(adjacency_matrix)
    # Create the network
    N = graph_from_data_frame(edge_list,directed=FALSE)
    m = ecount(N) # Number of edges
    n = vcount(N) # Number of nodes
    k = 2*m/n # Average degree
    apl = mean_distance(N,directed=FALSE) # Average path length
    c = transitivity(N) # Clustering coefficient of whole graph
    
    # Apply the function to the range of values and turn it into a dataframe
    df = data.frame(Rep,m,n,k,apl,c)
    colnames(df) = colnames(df.ER)
    df.ER=rbind(df.ER,df)
  }
  colnames(df.ER)=c("Rep","m","n","k","apl","c")
  return(df.ER)
}

ER = ER_function(k.ave,n.ave,edges,p)

#head(ER)
#head(df.reps)

# Check out average path length (See Connor et al. (2017))
plot(density(ER$apl),xlim=c(0,10),col="red")
lines(density(df.reps$apl))

# Check out clustering coefficient
plot(density(ER$c),xlim=c(0,0.6),col="red")
lines(density(df.reps$c))

# The red traces represent the distribution of values for a random network with similar properties
# We want our true values to lie outside of these; if the network is randomly all connected,
# average path length would be low because everything is connected. More structure increases APL.

# Since our network statistics seem to lie well outside of the null distributions 
# (more clustering and longer path length), it seems like there is a good chance that 
# our network is not just representing random noise. The clustering coefficient in the Connor paper was 0.38.

# Now we are ready to determine a consensus network. 
# Connor et al. ran 2000 network simulations with added noise and 
# using their permutationally-verified rho cutoff, and included edges 
# that were present in 90% of their simulations.
# We will do 1000 at 95%.
# Might recommend starting with fewer reps for computational tractability

#### Running consensus network ####

consensus_network = function(ps,rho){
  EL = data.frame(Rep=c(),X1=c(),X2=c())
  # Running it 1000 times 
  for (i in 1:1000){
    # Record iteration
    Rep = i
    # Making matrix of random values
    b = delta/1000
    E = replicate(dim(otu_table(ps))[1], rnorm(dim(otu_table(ps))[2]))
    E = 2*b*E
    E = -b + E
    # Add the noise to the matrix
    otu_table(ps) = otu_table(as.matrix(otu_table(ps))+abs(t(E)),taxa_are_rows=FALSE)
    # Calculate spearman correlations
    ps.dist = cor(otu_table(ps), use="everything", method="spearman") 
    # Turns Spearman into matrix
    adjacency_matrix = as.matrix(ps.dist)
    # Set correlations below rho cutoff to 0, and above to 1
    adjacency_matrix[abs(adjacency_matrix)<rho] = 0
    adjacency_matrix[abs(adjacency_matrix)>rho] = 1
    # Remove all taxa that have no correlations above the cutoff.
    adjacency_matrix = adjacency_matrix[rowSums(adjacency_matrix)>0,]
    # Create adjacency matrix from correlation matrix
    adjacency_matrix = graph.adjacency(adjacency_matrix, diag=FALSE, mode=c("undirected"), weighted=TRUE)
    # Gets edge list from adjacency matrix
    edge_list = get.edgelist(adjacency_matrix)
    edge_list = data.frame(Rep,edge_list)
    EL = rbind(EL,edge_list)
  }
  
  EL$X3=""
  EL$X3 = apply(EL,1,function(x) paste(sort(c(paste(x[2]),paste(x[3])))[1],sort(c(paste(x[2]),paste(x[3])))[2]))
  # The issue is that the pairing "OTU1 OTU2" can also be written "OTU2 OTU1"
  # If we want to match the pairings, we need to make sure they are always written the same way
  # To do this, we take X1 and X2 (the two nodes) and put them in a consistent order (~alphabetical)
  # (The above code is quite a bit (1-2x) faster than a paralellized mapply or a for loop to do this.)
  return(EL)
}

EL = consensus_network(ps,rho)

# Finds how many instances of each pairing there are
Common = EL %>%
  dplyr::group_by(X3)%>%
  dplyr::summarize(n())

# We want those that are present 95% of the time
# For 1000 runs, that's n>=950 or higher.
# For 100 runs, that's n>=95 or higher.
Consensus = Common[Common[,2]>=950,]

dim(EL)
dim(Common)
dim(Consensus)

# This is a good place to save!
# saveRDS(Consensus,"Derived_data/Network_CON18_4X")

#### Create iGraph of network, and identifying components of interest ####
## If you load a RDS/network from here, you also must create an appropriate ps object

# Divide the 95%-present paired responders back into two columns
# This is our consensus edge list
edge_list_consensus = separate(Consensus,1,into=c("X1","X2"),sep=" ")
edge_list_consensus = edge_list_consensus[,1:2]

# Setting edge list to the consensus value
edge_list = edge_list_consensus

colnames(edge_list)=c("A","B")
edge_list=data.frame(edge_list)
edge_list$CorVal=""
edge_list$CorSign=""

# Add correlation values and colors for the consensus edges
Cor_sign_estimate = function(ps,edge_list){
  ps.dist = cor(otu_table(ps), use="everything", method="spearman") 
  for (i in 1:dim(edge_list)[1]){
    CorVal = ps.dist[which(rownames(ps.dist)==data.frame(edge_list)$A[i]),
                     which(colnames(ps.dist)==data.frame(edge_list)$B[i])]
    CorSign = ifelse(CorVal>0,"Positive","Negative")
    edge_list[i,3]=CorVal
    edge_list[i,4]=CorSign
  }
  edge_list$EdgeColor[edge_list$CorSign=="Positive"] = "black"
  edge_list$EdgeColor[edge_list$CorSign=="Negative"] = "red"
  return(edge_list)
}

edge_list = Cor_sign_estimate(ps,edge_list)

igraph = graph_from_data_frame(edge_list,directed=FALSE)

# Get the taxon info for nodes
node_list = data.frame(V(igraph)$name)

# Add the OTU ID column
colnames(node_list)[1] = "OTU"

# Calculate standard properties of networks that are often reported in a typical study

network_properties = function(edge_list){
  N = graph_from_data_frame(edge_list,directed=FALSE)
  # Creates the network
  
  m = ecount(N) # Number of edges
  n = vcount(N) # Number of nodes
  k = 2*m/n # Average degree
  apl = mean_distance(N,directed=FALSE) # Average path length
  c = transitivity(N, type="global")
  cAve = transitivity(N, type = "average")
  # Clustering coefficient of whole graph - 
  # Transitivity measures the probability that the adjacent vertices of a vertex are connected.
  
  cl.mean = mean(closeness(N)) # closeness of graph (may not work if the graphs isn't fully connected)
  cl.sd = sd(closeness(N))
  ed = edge_density(N) # edge density of graph
  d = diameter(N) # diameter of graph

  # Turn it into a dataframe
  df = data.frame(m,n,k,apl,c,cAve,cl.mean,cl.sd,ed,d)
  return(df)
}

# Check out the full network properties. 
# These will change if/after you adjust for CorValue or remove unconnected components, below.
net.prop = network_properties(edge_list)
net.prop


# From https://chengjunwang.com/web_data_analysis/demo2_simulate_networks/
# plot and fit the power law distribution
fit_power_law = function(edge_list) {
  N = graph_from_data_frame(edge_list,directed=FALSE)
  # calculate degree
  d = degree(N, mode = "all")
  dd = degree.distribution(N, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
       col = 1, main = "Degree Distribution")
  curve(power.law.fit, col = "red", add = T, n = length(d))
}

fit_power_law(edge_list)


### Detect Modules
# Modules can be detected using cluster_spinglass
# Weighting by CorVal will positively weight co-occurences and negatively weight co-exclusions
Modules = cluster_spinglass(graph=igraph,weights=edge_list$CorVal,spins=100,implementation="neg")

# The cluster_spinglass will throw an error if there are unconnected components of the network
# If so, we can pare down graph so it's only got the connected vertices
# 2018 data--All networks have unconnected componenets, except for 32X

# Check out the various components (their sizes, membership)
components(igraph)

# Remove unconnected components from network
igraph2 = induced_subgraph(igraph, V(igraph)[components(igraph)$membership == which.max(components(igraph)$csize)])

# ## OR, select a specific component (comp==) to keep, if you want to analyze any component that isn't the biggest
# V(igraph)$comp = components(igraph)$membership
# igraph2 = induced_subgraph(igraph, V(igraph)$comp==1)

# Check the pre and post graphs
graphjs(igraph)
graphjs(igraph2)

# Update graph, edges and nodes (overwriting the full networks edge_list and node_list!)
igraph = igraph2
edge_list = edge_list[edge_list$A %in% names(V(igraph2)) & edge_list$B %in% names(V(igraph2)),]
node_list = data.frame(OTU=node_list[node_list$OTU %in% names(V(igraph2)),])

# Now, detect modules using cluster_spinglass
Modules = cluster_spinglass(graph=igraph,weights=edge_list$CorVal,spins=100,implementation="neg")

# Check out modules and modularity
Modules
modularity(Modules)

# Modularity (M) is an index measuring the extent to which a network is divided into modules, 
# and we used M > 0.4 as the threshold to define modular structures (Newman 2006)
# not all of our networks have meaningful modularity; some are <0.4. This makes sense for some treatments.

# Check out new network properties
net.prop = network_properties(edge_list)
net.prop

# Check out new power law distribution
fit_power_law(edge_list)

# Connectivity of each node can be determined based on its within-module connectivity (Zi)
# and among-module connectivity (Pi) (Guimera & Amaral 2005)

# Zi and Pi can be used to classify the nodes based on the topological roles they play in the network
# Node topologies are organised into four categories: module hubs (highly connected nodes within modules, Zi > 2.5)
# network hubs (highly connected nodes within entire network, Zi > 2.5 and Pi > 0.62)
# connectors (nodes that connect modules, Pi > 0.62)
# and peripherals (nodes connected in modules with few outside connections, Zi < 2.5 and Pi < 0.62) (Olesen et al. 2007; Zhou et al. 2010; Deng et al. 2012).

# To calculate the Zi and Pi of each node:
# First, find out which module it is in
# Make a list of all the other nodes in that module
# Calculate the connectivity of that node to all those other nodes
# Do this for each node
# Then, Zi is calculated as:
# (number of links from a given node to other nodes in the module - the average number for nodes in this module)
# Divided by the standard deviation of this value for nodes in this module.

# Then, do the same, but make the nodes list all the nodes in other modules

# Calculating Zi
adding_Zi = function(node_list,Modules,igraph){
  Zi=data.frame(Name=node_list$OTU,ModuleNumber=rep(0,length(node_list$OTU)),CON=rep(0,length(node_list$OTU)))
  for (i in 1:length(node_list$OTU)){
    node = paste(node_list$OTU[i])
    ModuleNumber = Position(function(x) node %in% x, Modules[], nomatch = 0)
    if(ModuleNumber !=0){
      ModuleNodes = Modules[[ModuleNumber]]
      modgraph = induced_subgraph(graph=igraph, v=c(ModuleNodes,node))
      CON=try(as.numeric(lengths(adjacent_vertices(modgraph,node))),TRUE)
      if(isTRUE(class(CON)=="try-error")) { CON=NA } else {CON = as.numeric(lengths(adjacent_vertices(modgraph,node)))}
    }
    if(ModuleNumber ==0){CON=NA}
    Zi$Name[i]=as.factor(node_list$OTU[i])
    Zi$ModuleNumber[i]=ModuleNumber
    Zi$CON[i]=CON
  }
  
  # (number of links from a given node to other nodes in the module - the average number for nodes in this module)
  # Divided by the standard deviation of this value for nodes in this module.
  Zi = Zi %>%
    filter(!is.na(CON))%>%
    group_by(ModuleNumber)%>%
    mutate(MeanCON=mean(CON))%>%
    mutate(SdCON=sd(CON))%>%
    mutate(Zi=((CON-MeanCON)/SdCON))
  return(Zi)
}
Zi = adding_Zi(node_list,Modules,igraph)

# Next, we add Pi
adding_Pi = function(node_list,Modules,igraph){
  Pi=data.frame(Name=rep(node_list$OTU,dim(Modules[])),HomeModuleNumber=rep(0,length(node_list$OTU)),OtherModuleNumber=rep(0,length(node_list$OTU)),TotalCON=rep(0,length(node_list$OTU)),CON=rep(0,length(node_list$OTU)))
  for (i in 1:length(node_list$OTU)){
    node = paste(node_list$OTU[i])
    HomeModuleNumber = Position(function(x) node %in% x, Modules[], nomatch = 0)
    ModuleNumbers = 1:dim(Modules[])
    n = length(ModuleNumbers)
    TotalCON = as.numeric(lengths(adjacent_vertices(igraph,node)))
    lowend = (i-1)*n+1
    highend = n*i
    Pi$Name[lowend:highend]=node_list$OTU[i]
    Pi$HomeModuleNumber[lowend:highend]=HomeModuleNumber
    Pi$OtherModuleNumber[lowend:highend]=ModuleNumbers
    Pi$TotalCON[lowend:highend]=TotalCON
    if(HomeModuleNumber !=0){    
      for (j in ModuleNumbers){
        OtherModuleNumber = j
        NodesInOtherModule = Modules[[OtherModuleNumber]]
        modgraph = induced_subgraph(graph=igraph, v=c(node,NodesInOtherModule))
        CON=as.numeric(lengths(adjacent_vertices(modgraph,node)))
        Pi$CON[Pi$HomeModuleNumber==HomeModuleNumber & Pi$OtherModuleNumber==OtherModuleNumber & Pi$Name==node]=CON
      }
    }
  }
  Pi$kk2 = (Pi$CON/Pi$TotalCON)^2
  return(Pi)
}

Pi = adding_Pi(node_list,Modules,igraph)

Pifinal = Pi %>%
  dplyr::group_by(Name,HomeModuleNumber,TotalCON)%>%
  dplyr::summarize(Sum=sum(kk2))%>%
  dplyr::mutate(Pi=1-Sum)

# Calculating Zi from Pi data
Zinew = Pi %>% filter(HomeModuleNumber==OtherModuleNumber) %>% mutate(MeanCON=mean(CON),SdCON=sd(CON),Zi=((CON-MeanCON)/SdCON))

# Bringing module data together
# including identification of hubs and connectors within the network.
# Thresholds based on Modules paper
Making_module_data = function(Pifinal,Zinew){
  Pthresh = 0.62
  Zthresh = 2.5
  ModuleData=data.frame(Name=Pifinal$Name,Module=Pifinal$HomeModuleNumber,TotalCON=Pifinal$TotalCON,ModuleCON=Zinew$MeanCON,Pi=Pifinal$Pi,Zi=Zinew$Zi)
  ModuleData$Class = ifelse(ModuleData$Zi>Zthresh & ModuleData$Pi>Pthresh,"Network Hub",
                            ifelse(ModuleData$Zi>Zthresh & ModuleData$Pi<Pthresh,"Module Hub",
                                   ifelse(ModuleData$Zi<Zthresh & ModuleData$Pi>Pthresh,"Connector", "Peripheral")))
  return(ModuleData)
}

ModuleData = Making_module_data(Pifinal,Zinew)

p = ggplot(ModuleData)
p = p + geom_point(aes(x=Pi,y=Zi,color=Class))
p
tail(ModuleData)


## Add this info to the nodes list
add_modInfo = function(node_list,ModuleData){
  node_list$Pi=c()
  for (i in 1:dim(node_list)[1]){
    OTU = node_list$OTU[i]
    x = ifelse(OTU %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==OTU,]$Pi!=(-Inf),ModuleData[ModuleData$Name==OTU,]$Pi, NA))
    node_list$Pi[i] = x
  }
  
  node_list$Zi=c()
  for (i in 1:dim(node_list)[1]){
    OTU = node_list$OTU[i]
    x = ifelse(OTU %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==OTU,]$Zi!=(-Inf),ModuleData[ModuleData$Name==OTU,]$Zi, NA))
    node_list$Zi[i] = x
  }
  
  node_list$NetworkRole=c()
  for (i in 1:dim(node_list)[1]){
    OTU = node_list$OTU[i]
    x = ifelse(OTU %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==OTU,]$Class!=(-Inf),ModuleData[ModuleData$Name==OTU,]$Class, NA))
    node_list$NetworkRole[i] = x
  }
  
  node_list$Module=c()
  for (i in 1:dim(node_list)[1]){
    OTU = node_list$OTU[i]
    x = ifelse(OTU %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==OTU,]$Module!=(-Inf),ModuleData[ModuleData$Name==OTU,]$Module, NA))
    node_list$Module[i] = x
  }

    node_list$NetworkRoleColour = ifelse(node_list$NetworkRole=="Connector","red",
                                       ifelse(node_list$NetworkRole=="Module Hub","navy","white"))
  
  return(node_list)
}

node_list = add_modInfo(node_list,ModuleData)


# Check most abundant modules
ModProps = node_list %>% dplyr::group_by(Module) %>% dplyr::summarize(Total=n()) %>% dplyr::arrange(-Total)
head(ModProps)
dim(ModProps)

#### Adding taxonomy to the nodes list, in order to readily see which nodes were which taxa ####

# Fix up some taxon naming issues
ignoreList = c("","uncultured","uncultured bacterium","uncultured soil bacterium","uncultured forest soil bacterium", "Ambiguous_taxa",
               "uncultured actinobacterium","uncultured planctomycete","uncultured Chloroflexi bacterium","uncultured Syntrophobacteraceae bacterium",
               "uncultured crenarchaeote","uncultured Gemmatimonadetes bacterium", "uncultured Acidobacteria bacterium",
               "uncultured Planctomyces sp.", "metagenome", "Subgroup 6", "Blastocatellia (Subgroup 4)", "uncultured Holophaga sp.",
               "uncultured Hyphomicrobiaceae bacterium", "uncultured proteobacterium")


taxonomy_labeller = function(ps.trim,node_list){
  node_list$Species = c()
  node_list$Genus = c()
  node_list$Family = c()
  node_list$Order = c()
  node_list$Class = c()
  node_list$Phylum = c()
  # node_list$Kingdom = c()
  for (i in 1:dim(node_list)[1]){
    OTU = node_list$OTU[i]
    Species = ifelse(OTU %in% row.names(tax_table(ps.trim)),paste(tax_table(ps.trim)[paste(OTU),"Species"]),"")
    Genus = ifelse(OTU %in% row.names(tax_table(ps.trim)),paste(tax_table(ps.trim)[paste(OTU),"Genus"]),"")
    Family = ifelse(OTU %in% row.names(tax_table(ps.trim)),paste(tax_table(ps.trim)[paste(OTU),"Family"]),"")
    Order = ifelse(OTU %in% row.names(tax_table(ps.trim)),paste(tax_table(ps.trim)[paste(OTU),"Order"]),"")
    Class = ifelse(OTU %in% row.names(tax_table(ps.trim)),paste(tax_table(ps.trim)[paste(OTU),"Class"]),"")
    Phylum = ifelse(OTU %in% row.names(tax_table(ps.trim)),paste(tax_table(ps.trim)[paste(OTU),"Phylum"]),"")
    node_list$Species[i] = Species
    node_list$Genus[i] = Genus
    node_list$Family[i] = Family
    node_list$Order[i] = Order
    node_list$Class[i] = Class
    node_list$Phylum[i] = Phylum
    # node_list$Kingdom[i] = Kingdom
    node_list$SpecialName[i] = paste(ifelse(Genus %in% ignoreList | is.na(Genus),
                                         ifelse(Family %in% ignoreList | is.na(Family),
                                                ifelse(Class %in% ignoreList | is.na(Class),
                                                       ifelse(Phylum %in% ignoreList | is.na(Phylum),paste(OTU),
                                                              paste(Phylum)),paste(Class)),paste(Family)),paste(Genus)))
  }
  return(node_list)
}
node_list = taxonomy_labeller(ps.trim,node_list)


##### Making network figures ######

set.seed(100)

# Edge colours
E(igraph)$color = c(edge_list$EdgeColor)
E(igraph)$size = 4

# Plot the 3D network graph
graphjs(igraph,
        width=700,height=500,vertex.shape="circle",bg="white",fg="white",vertex.size=0.5,
        edge.alpha=0.75)


#### Modules figure ####
set.seed(4)

igraph.modules = list()
for (i in 1:length(levels(as.factor(node_list$Module)))){
  Mod = levels(as.factor(node_list$Module))[i]
  igraph.modules[[i]] = induced_subgraph(graph=igraph, v=c(Modules[[Mod]]))
}

# get the coordinates of each individual module graph
graphs = igraph.modules
# Put each module's nodes into a circle
layouts = lapply(graphs,layout_in_circle)
# take our module graphs, and their layouts, and put them all in the same space
# They are placed in order of size
# This results in a long matrix with the position of each vertex
lay = merge_coords(graphs, layouts)
# In this step, the graphs are all merged
# It assumes no overlap (or relabels vertices to create this)
g = disjoint_union(graphs)

# Take our set of coordinates
layOrder = data.frame(lay)
# Then we assume they are in the same/correct order as in the graph list
layOrder$OTU = V(g)$name
# Here, we then put them in the same order as the igraph names
layOrder$OTU = factor(layOrder$OTU, levels = V(igraph)$name)

layOrder = layOrder%>%
  dplyr::arrange(OTU)
lay = as.matrix(layOrder[,1:2])

# Go back to the original network, and give it this layout
igraph.Modules = igraph
igraph.Modules$layout = lay
# This tells us where each node should be

plot.igraph(igraph.Modules, vertex.size=3, vertex.label=NA,
            edge.width=c(0.2),edge.alpha=c(0.1))


#### If we wanted to select certain modules... ####

### Pull out all modules with more than (modcutoff) nodes
# modcutoff = 3
# SelectMods = ModProps %>% filter(Total>modcutoff) %>% select(Module)
# SelectMods = SelectMods$Module
# SelectMods
### OR: Identify certain modules of interest
SelectMods = c(1,2,3)
set.seed(4)

#Pull out selected modules
igraph.3Modules = induced_subgraph(igraph, node_list$OTU[node_list$Module %in% c(SelectMods)])
node_list.3Modules = node_list[node_list$Module %in% c(SelectMods),]

# Make list for only the modules of interest
igraph.modules.3Modules = list()
for (i in 1:length(levels(as.factor(node_list.3Modules$Module)))){
  Mod = levels(as.factor(node_list.3Modules$Module))[i]
  igraph.modules.3Modules[[i]] = induced_subgraph(graph=igraph.3Modules, v=c(Modules[[Mod]]))
}
# Same code as above for plotting
graphs = igraph.modules.3Modules
layouts = lapply(graphs,layout_in_circle)
lay = merge_coords(graphs, layouts)
g = disjoint_union(graphs)
layOrder = data.frame(lay)
layOrder$OTU = V(g)$name
layOrder$OTU = factor(layOrder$OTU, levels = V(igraph.3Modules)$name)
layOrder = layOrder%>%
  dplyr::arrange(OTU)
lay = as.matrix(layOrder[,1:2])
igraph.3Modules$layout = lay
plot.igraph(igraph.3Modules, vertex.size=3, vertex.label=NA,
            edge.width=c(0.2),edge.alpha=c(0.1))




#### Cross-module properties by mixing treatment ####

ps.TRT.norm2 = subset_samples(ps.TRT.norm, VortexControl == "N")
ps.TRT.norm2 = prune_taxa((taxa_sums(ps.TRT.norm2) > 0), ps.TRT.norm2)

mdf = psmelt(ps.TRT.norm2)

mdf2 = mdf %>%
  dplyr::filter(OTU %in% node_list$OTU)

# Get only network OTUs, change naming scheme to match node_list
M = node_list[,c("OTU","Module")]
# Get just the two relevant columns
mdf2 = merge(mdf2,M,by="OTU",all.x=TRUE)
# Merge the two dataframes

OTUbyMod = mdf2 %>%
  dplyr::group_by(Sample,Module,TimesMixed,)%>%
  dplyr::summarize(Abundance = sum(Abundance))%>%
  dplyr::mutate(ModuleName = paste("Module",Module))

## Option to choose which modules to include
#plotmods = c(1,2)
#plotmods = SelectMods

d = OTUbyMod %>%
  #dplyr::filter(Module %in% plotmods)%>%
  dplyr::mutate(StripName = paste(Module))%>%
  dplyr::group_by(Sample,Module,TimesMixed)%>%
  dplyr::summarize(Abundance=sum(Abundance))

ModLevels = ModProps$Module # for putting modules in order of most OTUs to fewest
d=d%>%
  mutate(highlight=ifelse(TimesMixed=="4","Highlighted","Normal"))

#modtheMod=as_labeller(c('1'='6', '2'='7', "4" = "8", "3" = "9")) #to modify module numbers in facet labels; must include labeller=modtheMod in facet_wrap()

p = ggplot(data=d, aes(x=TimesMixed,y=Abundance, fill=highlight))
p = p + geom_boxplot() +
  scale_fill_manual(values=c("#69b3a2","white"))
p = p + theme_bw()
p = p + facet_wrap(~factor(Module, levels=ModLevels), scales="free_y",ncol=6) + expand_limits(y=0)
p = p + ylab("Relative Abundance") + xlab("Mixing frequency")
p = p + theme(axis.text=element_text(size=10),
              axis.title=element_text(size=12),
              strip.text=element_text(size=12),
              axis.text.x=element_text(angle=45, vjust=1, hjust=1))+
  theme(legend.position = "None")
p 

#### Cross-module properties, WITHIN mixing treatment, by mixing set ####
#mdf = psmelt(ps.TRT.norm)
mdf3 = mdf %>%
  dplyr::filter(OTU %in% node_list$OTU)

# Get only network OTUs, change naming scheme to match node_list
M = node_list[,c("OTU","Module")]
# Get just the two relevant columns
mdf3 = merge(mdf3,M,by="OTU",all.x=TRUE)
# Merge the two dataframes

mdf3 = subset(mdf3, VortexControl=="N" & TimesMixed=="4")

OTUbyMod3 = mdf3 %>%
  dplyr::group_by(Sample,Module,Rep, TimesMixed)%>%
  dplyr::summarize(Abundance = sum(Abundance))%>%
  dplyr::mutate(ModuleName = paste("Module",Module))

# Option to choose which modules to include 
# plotmods = c(1,2)
# plotmods = SelectMods

d3 = OTUbyMod3 %>%
  #dplyr::filter(Module %in% plotmods)%>%
  dplyr::mutate(StripName = paste(Module))%>%
  dplyr::group_by(Sample,Module,Rep,TimesMixed)%>%
  dplyr::summarize(Abundance=sum(Abundance))
ModLevels = ModProps$Module

p = ggplot(data=d3, aes(x=Rep,y=Abundance,color=Rep))
#p = p + geom_jitter(alpha = 0.4) # adds scatter
p = p + geom_boxplot(alpha = 0.2)
p = p + theme_bw()
p = p + scale_color_manual(values=wes_palette(4, name = "GrandBudapest1", type = "continuous"), name = "")
p = p + facet_wrap(~factor(Module, levels=ModLevels), labeller=modtheMod, scales="free_y",ncol=6) + expand_limits(y=0)
p = p + ylab("Relative Abundance") + xlab("Mixing Set")# + guides(color = "none")
p = p + theme(axis.text=element_text(size=10),
              axis.title=element_text(size=12),
              strip.text=element_text(size=12))
p = p + theme(legend.position = "None")
p 


# # Plot for 1x, Initial (no reps/mixing sets)
# p = ggplot(data=d3, aes(x=TimesMixed, y=Abundance))
# p = p + geom_jitter(alpha = 0.4)
# p = p + geom_boxplot(alpha = 0.2)
# p = p + theme_bw()
# p = p + scale_color_manual(values=wes_palette(4, name = "GrandBudapest1", type = "continuous"), name = "")
# p = p + facet_wrap(~factor(Module, levels=ModLevels), scales="free_y",ncol=5) + expand_limits(y=0)
# p = p + ylab("Relative Abundance") + xlab("")
# p = p + theme(axis.text=element_text(size=10),
#               axis.title=element_text(size=12),
#               strip.text=element_text(size=12))
# p = p + theme(legend.position = "None")
# p #750 x 500

