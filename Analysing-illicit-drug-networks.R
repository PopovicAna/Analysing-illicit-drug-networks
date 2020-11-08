source("Files/Prioritisation-of-analytical-techniques.R")

library(lubridate)
library(igraph)
library(DescTools)

# Snippet of code from Dual-apporach-for-score evaluation -----------------

# Defining TP, TN, FP and FN rates at each threshold values (THV - aka comparison metric score)
Cutoffs <-  data.frame(
  THV = RCM(OPT_GCMS_CM_R)$thresholds, 
  TPR = RCM(OPT_GCMS_CM_R)$sensitivities, 
  FPR = 1-RCM(OPT_GCMS_CM_R)$specificities,
  FNR = 1-RCM(OPT_GCMS_CM_R)$sensitivities, 
  TNR = RCM(OPT_GCMS_CM_R)$specificities)

# Setting an acceptable FPR to define whether specimen pairs are linked or not
# Note: For the purpose of this research an acceptable FPR is 0.025 (i.e. 2.5%)
LINK_THV <- as.numeric(tail(subset(Cutoffs, FPR > 0.025, select = THV), n = 1))



# Creating clusters or Chemical Classes (CC) for specimens ----------------

# Calculating similarity between specimens based on the optimal comparison metric (CM)
Scores_GCMS <- as.dist(get(OPT_GCMS_CM)(as.matrix(get(names(OPT_AT)))), diag = FALSE, upper = FALSE)

# Defining the Optimum linkage method for clustering analysis
OPT_LM <- rownames(as.matrix(which.max(mapply(
  function(LM){cor(Scores_GCMS, cophenetic(hclust(d = Scores_GCMS, method = LM)))},
  c("ward.D","ward.D2","single","complete","average","mcquitty")))))

# Defining chemical classes using hierarchical clustering analysis and the THV (cutoff value for clustering)
CC_GCMS <- hclust(Scores_GCMS, method = OPT_LM)
clusters <- data.frame(CC = dendextend::cutree(CC_GCMS, h = LINK_THV, order_clusters_as_data = F))

"Include code from SimAPP for 2PGs?"

# Adding specimen-CC codes to Lookup
Lookup$cluster <- clusters$CC[match(Lookup$Specimen, rownames(clusters))]

# Visualising dendrogram of chemical classes
plot(CC_GCMS, xlab = "Specimens", ylab = "CM Score", sub = "", labels = FALSE, hang = -1)
abline(h = LINK_THV, lty = 2)

# Visualising an example of one Chemical class
Sub_Dend=function(n){plot(cut(as.dendrogram(CC_GCMS),h=LINK_THV)$lower[[n]],
                          main = paste0("CC (",n,") of main dendrogram with THV at h=",round(LINK_THV,2)),
                          xlab = "Specimens", ylab = "CM Score", )
} 
Sub_Dend(6)


# Creating network plots between CCs and respective Groups ----------------
#Note: Trying to plot the relationship betweeen CC and specimen results in a messy graph.
#A cleaner alternative is plotting the relationships between CC and specimen groups. 

# Extracting date information (onset and terminus) for each CC
Net_CC <- do.call(data.frame, aggregate(Date~cluster,Lookup, function(x) c(min(x),max(x))))
Net_CC <- setNames(Net_CC,c("vertex.id","onset","terminus"))
Net_CC$terminus <- Net_CC$terminus+1
Net_CC$type <- "CC"
Net_CC$name <- paste0("CC_",Net_CC$vertex.id)

# Extracting date information (onset and terminus) for each specimen group (SG)
Net_SG <- Lookup[,c("Group","cluster")]
names(Net_SG)[1] <- "name"
Net_SG$terminus <- Net_CC$terminus[match(Net_SG$cluster,Net_CC$vertex.id)]
Net_SG <- aggregate(terminus~name,Net_SG,max)
Net_SG$onset <- as.numeric(Lookup$Date[match(Net_SG$name,Lookup$Group)])
Net_SG$type <- "SG"
Net_SG$vertex.id <- Net_SG$name

# Creating a dataframe of unique nodes
Nodes <- rbind(Net_CC,Net_SG)
Nodes$onset.censored <- F
Nodes$terminus.censored <- F
Nodes <- Nodes[c("name","onset","terminus","vertex.id","onset.censored","terminus.censored","type")]

# Creating a dataframe of unique links
Links <- data.frame(from=Lookup$cluster,to=Lookup$Group, stringsAsFactors = F)
Links$onset <- Nodes$onset[match(Links$to,Nodes$vertex.id)]
Links$terminus <- Nodes$terminus[match(Links$from,Nodes$vertex.id)]
Links$onset.censored <- F
Links$terminus.censored <- F
Links$from <- paste0("CC_",Links$from)

# Creating an undirected graph 
Net <- graph_from_data_frame(d = Links[,1:5],vertices = Nodes[,1:3], directed = F)
plot(Net)

# Plotting CC-Group networks between specific dates (e.g. first year of data)
Sub_Net <- delete_edges(Net, which(E(Net)$onset<as.numeric(as.Date("2016-01-01")) |   
                                     E(Net)$terminus>as.numeric(as.Date("2016-12-31"))))
Sub_Net <- delete_vertices(Sub_Net,degree(Sub_Net)==0)
plot(Sub_Net)



# Relational analysis: ----------------------------------------------------

# Number of specimens, groups and CCs
No_S <- length(unique(Lookup$Specimen))
No_SG <- length(unique(Lookup$Group))
No_CC <- length(unique(Lookup$cluster))

# Network density (bipartite network)
Net_density <- ecount(Net)/(
  length(names(V(Net))[substr(names(V(Net)),1,2)=="CC"])*
    (length(names(V(Net)))-length(names(V(Net))[substr(names(V(Net)),1,2)=="CC"]))
)

# CC connected to 2 or more Groups
CC_2_SG <- Lookup %>% 
  group_by(cluster) %>% 
  summarise(G = n_distinct(Group)) %>% 
  filter(G >= 2) %>% 
  count()

# Groups connected to 1  CC (Mono-profile Groups)
MOP <- Lookup %>% 
  group_by(Group) %>% 
  summarise(CC = n_distinct(cluster)) %>% 
  filter(CC == 1) %>% 
  count()
MOP_per <- (MOP/No_SG)*100

# Groups connected to 2 or more CCs (Multi-profile Groups)
MLP <- Lookup %>% 
  group_by(Group) %>% 
  summarise(CC = n_distinct(cluster)) %>% 
  filter(CC >= 2) %>% 
  count()
MLP_per <- (MLP/No_SG)*100

# Extracting the largest component in the network 
Net_Comps <- components(Net)
Net_LC <- decompose(Net)[[which.max(Net_Comps$csize)]]
plot(Net_LC)



# Temporal analysis: ------------------------------------------------------

# Number of CC gained, retained and lost per quarter
CC_Life <- Lookup %>% group_by(cluster) %>% summarise(min = min(Date), max = max(Date))
CC_Life$days <- difftime(CC_Life$max ,CC_Life$min , units = c("days")) + 1
CC_Life$month <- round(CC_Life$days/30,0)
CC_Life$int <- interval(CC_Life$min, CC_Life$max)

Quarters <- data.frame(Start = seq(as.Date("2015-01-01"), as.Date("2019-10-01"), by = "quarter"),
                       End = seq(as.Date("2015-04-01"), as.Date("2020-01-10"), by = "quarter")-1)
Quarters$int <- interval(Quarters$Start, Quarters$End)
Quarters$CC <- unlist(lapply(c(1:nrow(Quarters)),function(x){length(na.omit(intersect(CC_Life$int, Quarters$int[x])))}))
Quarters$gain <- unlist(lapply(c(1:nrow(Quarters)),function(x){sum(as.Date(CC_Life$min) %within% Quarters$int[x], na.rm = TRUE)}))
Quarters$const <- Quarters$CC-Quarters$gain
Quarters$lost <- unlist(lapply(c(1:nrow(Quarters)),function(x){sum(as.Date(CC_Life$max) %within% Quarters$int[x-1], na.rm = TRUE)}))
Quarters$Q <- quarters(Quarters$Start)
Quarters <- pivot_longer(Quarters, gain:lost)

# Visualising the number of CC gained, retained and lost per quarter 
Quarters %>% 
  mutate(
    value = case_when(
    name == "lost" ~ -1 * value,
    TRUE ~ as.numeric(value)),
    YQ = paste0(year(Start),"-",Q)
    ) %>% 
  ggplot(aes(YQ, value)) +
  geom_bar(aes(fill = name), stat = "identity") +
  scale_fill_manual(values = c("grey60","coral2", "cadetblue")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 90)
  )

# OTHER TEMPORAL
#mean(CC_Life$days)
#median(CC_Life$days)