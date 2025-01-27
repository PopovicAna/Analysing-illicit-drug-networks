source("Files/Prioritisation-of-analytical-techniques.R")

library(lubridate)
library(igraph)
library(DescTools)
library(plotly)

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
Sub_Dend(2)


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

# Defining node positions and colours for visualisation
LN <- as.data.frame(
  layout_with_fr(Net,
                 grid="nogrid",
                 coords=norm_coords(layout_with_fr(Net, grid = "nogrid"),
                                    xmin = -1, xmax = 1, ymin = -1, ymax = 1),  
                 niter=10,start.temp=0.05),
  row.names = names(V(Net))
)
LN$colour <-  ifelse(substr(rownames(LN),1,2)=="CC","black","grey")

# Extracting the network edgelist
es <- as.data.frame(get.edgelist(Net))

# Creating the network plot
plot_ly(
  x = LN[,1], y = LN[,2], type = "scatter", mode = "markers", source="subset",
  name = rownames(LN), text = rownames(LN), hoverinfo = rownames(LN),  
  marker = list(color = LN$colour)
) %>%
  layout(
    shapes = lapply(c(1:length(es[1]$V1)),function(i){
      list(
        type = "line",
        line = list(color = "#030303", width = 0.3),
        x0 = LN[rownames(LN)==as.character(es[i,]$V1),1],
        y0 = LN[rownames(LN)==as.character(es[i,]$V1),2],
        x1 = LN[rownames(LN)==as.character(es[i,]$V2),1],
        y1 = LN[rownames(LN)==as.character(es[i,]$V2),2]
      )
    }),
    xaxis = list(title = "", showgrid = F, showticklabels = F, zeroline = F),
    yaxis = list(title = "", showgrid = F, showticklabels = F, zeroline = F))


# Plotting CC-Group networks between specific dates (e.g. first year of data)
Sub_Net <- delete_edges(Net, which(E(Net)$onset<as.numeric(as.Date("2015-01-01")) |   
                                     E(Net)$terminus>as.numeric(as.Date("2015-12-31"))))
Sub_Net <- delete_vertices(Sub_Net,degree(Sub_Net)==0)
plot(Sub_Net)



# Relational analysis: ----------------------------------------------------

# Number of seizures and CCs
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
Net_Comps <- components(Net)$csize
Net_Comps_Order <- order(Net_Comps, decreasing = T)
Net_Dec <- decompose(Net)
Net_LC <- decompose(Net)[[which.max(Net_Comps)]]

# Extracting the 6 largest components in the network
Y <- lapply(1:6, function(x)Net_Dec[[Net_Comps_Order[x]]])

# Defining node positions and colours of the lergest components for visualisation
LN_Sub <- lapply(1:6,function(x)
  data.frame(layout_with_fr(Y[[x]]),
             colour = LN[match(names(V(Y[[x]])),rownames(LN)),3],
             row.names = names(V(Y[[x]]))))

# Extracting the edgelists of the largest components
es_Sub <- lapply(1:6,function(x)as.data.frame(get.edgelist(Y[[x]])))

# Visualising the largest components
LC_Plots <- lapply(
  1:6,
  function(x){
    plot_ly(
      x = LN_Sub[[x]][,1], y = LN_Sub[[x]][,2], type = "scatter", mode = "markers",
      name = rownames(LN_Sub[[x]]), text = rownames(LN_Sub[[x]]), hoverinfo = rownames(LN_Sub[[x]]),
      marker=list(color = LN_Sub[[x]]$colour
      )) %>%
      layout(
        plot_bgcolor='#FAFBFD',
        shapes = lapply(c(1:length(es_Sub[[x]][1]$V1)),function(i){
          list(
            type = "line",
            line = list(color = "#030303", width = 0.3),
            x0 = LN_Sub[[x]][rownames(LN_Sub[[x]])==as.character(es_Sub[[x]][i,]$V1),1],
            y0 = LN_Sub[[x]][rownames(LN_Sub[[x]])==as.character(es_Sub[[x]][i,]$V1),2],
            x1 = LN_Sub[[x]][rownames(LN_Sub[[x]])==as.character(es_Sub[[x]][i,]$V2),1],
            y1 = LN_Sub[[x]][rownames(LN_Sub[[x]])==as.character(es_Sub[[x]][i,]$V2),2]
          )
        }),
        xaxis = list(title = "", showgrid = F, showticklabels = F, zeroline = F),
        yaxis = list(title = "", showgrid = F, showticklabels = F, zeroline = F))
  })

subplot(LC_Plots, nrows = as.integer(sqrt(6)))




# Temporal analysis: ------------------------------------------------------

# Creating a line plot (with a cumulative reference) showing the number of months CCs are active
# Calculating the number of CCs observed over time
A <- aggregate(Date~cluster, Lookup, function(x){c(min(x), max(x))})
A <- as.data.frame(do.call(cbind,A))
A$V2 <- as.Date(A$V2, origin = "1970-01-01")
A$V3 <- as.Date(A$V3, origin = "1970-01-01")
A$days <- difftime(A$V3 ,A$V2 , units = c("days")) + 1
A$days <- as.numeric(A$days)
A$mo <- round(A$days/30,0)
A <<- A

B <- aggregate(cluster~mo, A, function(x){length(unique(x))})
B$cumu <- cumsum(B$cluster)
B$cumu <- round(B$cumu/max(B$cumu)*100,0)

# visualising the no of CC observed over time 
ggplot(B,aes(x = mo)) + 
  geom_line(aes(y=cluster),size=1.5)+
  geom_line(aes(y=cumu/2),size=1.5, color="grey60")+
  scale_x_continuous(name="Number of months",expand = c(0,0), limits = c(0,60))+
  scale_y_continuous(name="Number of CCs",breaks=seq(0,50,5),
                     sec.axis=sec_axis(~.*2,name="Cumulative number of CCs (%)",breaks=seq(0,100,20)))+
  theme_light()+
  theme(legend.position="bottom",
        axis.title=element_text(face="bold", size = 12), 
        axis.text = element_text(size = 10),
        strip.text=element_text(face="bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

# mean and median number of days a CC has been observed 
round(mean(A$days),0)
round(median(A$days),0)



# Plotting Cc 'life' over all quarters 
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
    axis.text.x = element_text(angle = 90, size = 10),
    axis.text.y = element_text(size = 10)
  )



# Spatial analysis: ------------------------------------------------------
library(rgdal)
library(leaflet)
library(htmltools)
# Importing Australian postcode shape data
PC <- readOGR(dsn = "Data/PC_Shapes/V2/POA_2016_AUST_V2.shp")

# creating a map to display where seizures have been made (respective to the date slider)
# NOTE: CURRENTLY THIS ONLY EXTENDS TO NSW DAATA
countPC <-
  Lookup %>% 
  filter(!is.na(Postcode)) %>%
  filter(!is.na(cluster)) %>%
  select(Group,Postcode,Date) %>%
  distinct()

# Filtering seizure data for certain dates
data_input <-
  countPC %>%
  filter(Date >= "2015-01-01"
  ) %>%
  filter(Date <= "2019-12-31"
  ) %>%
  group_by(Postcode) %>%
  summarise(Seizures = dplyr::n())

# Extracting the subset of postcodes respective to the chosen date range
PCs <- subset(PC, is.element(PC$POA_NAME16,data_input$Postcode))

# ordering the data based on postcodes
data_input_ordered <- 
  data_input[order(match(data_input$Postcode,PCs$POA_NAME16)),]

# creating labels for teh postcode polygons
labels <- paste("<p>", "Postcode: ", data_input_ordered$Postcode, "<p>",
                "<p>", "No. of seizures: ", data_input_ordered$Seizures, "<p>", sep = "")

# Creating a colour palatte for number of seizures
pal <- colorBin("YlOrRd", domain = c(0,1), bins = seq(0,8,2))

# Plotting the seizures per postode map
leaflet() %>% 
  setView(lat = -33, lng = 147, zoom = 6) %>%
  addProviderTiles(providers$Stamen.TonerLite) %>%    # Map Themes
  addPolygons(data = PCs, 
              color = "#666666", 
              weight = 1, 
              fillOpacity = 0.8, 
              fillColor = pal(data_input_ordered$Seizures),
              highlightOptions = highlightOptions(
                weight = 5,
                color = "#666666",
                fillOpacity = 0.7,
                bringToFront = TRUE),
              label = lapply(labels, HTML)) %>%
  addLegend(title = "No. of seizures",
            pal = pal,
            values = data_input_ordered$Seizures, 
            opacity = 0.7,
            position = "topright")



# Quantitative analysis: ------------------------------------------------------

# Calculating purities for each Region and all regions (domestic) per quarter
  Lkp_Regional <- Lookup
  Lkp_Domestic <- Lookup
  Lkp_Domestic$Region = "Domestic"
  
  Lkp_All <-  rbind(Lkp_Regional[,c("Specimen","Date","Year","Purity","Pre","Region")],
                    Lkp_Domestic[,c("Specimen","Date","Year","Purity","Pre","Region")])
  Lkp_All$Quarter <- quarters(Lkp_All$Date)
  Lkp_All$Region_f <- factor(Lkp_All$Region, 
                             levels = c("ACT","NSW","NT","QLD","SA","TAS","VIC","WA","Domestic"))
  
  ggplot(Lkp_All, aes(x=Quarter, y=Purity)) +
    geom_boxplot(outlier.shape = 21) + 
    facet_grid(Region_f~Year, scales="free_x", drop=T) +
    labs(y="Purity (%)")+
    theme_light()+
    theme(axis.title=element_text(face="bold",size = 11),
          axis.text = element_text(size = 10),
          strip.text=element_text(face="bold",size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
# Determining precursors used for each Region and all regions (domestic)
  Lkp_All$Pre <- ifelse(is.na(Lkp_All$Pre),"Unclassified",
                        ifelse(Lkp_All$Pre%in%c("Pred. EPH/PSE","PSE","EPH"),"EPH/PSE",
                               ifelse(Lkp_All$Pre=="Pred. P2P","P2P",Lkp_All$Pre)))
  
  Lkp_Pre <- Lkp_All %>% 
    group_by(Region,Year) %>% 
    count(Pre) %>% 
    group_by(Region,Year) %>% 
    mutate(countT = sum(n)) %>% 
    group_by(Pre, add=TRUE) %>%
    mutate(per=round(100*n/countT,2))
  Lkp_Pre$Region_f <- factor(Lkp_Pre$Region, 
                             levels = c("ACT","NSW","NT","SA","TAS","VIC","WA","Domestic","National"))
  Lkp_Pre$Pre_f <- factor(Lkp_Pre$Pre, 
                          levels = c("EPH/PSE","P2P","Mixed","Unclassified"))
  
  ggplot(Lkp_Pre,aes(x=Pre_f,y=per)) +
    geom_bar(stat="identity")+
    facet_grid(Region_f~Year, scales="fixed", drop=T) +
    labs(y="Percentage", x="Precursor Type")+
    scale_fill_grey()+
    theme_light()+
    theme(axis.title=element_text(face="bold",size = 11),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(angle = 45,hjust=1),
          strip.text=element_text(face="bold",size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
