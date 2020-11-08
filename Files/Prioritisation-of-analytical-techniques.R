
# Reading in the data and basic tidying: ----------------------------------

# Load the required libraries
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(ggcorrplot)
library(pROC)
library(ggh4x)

# Importing the lookup table
Lookup <-  read.xlsx("Data/Profiles.xlsx", sheet = "Lookup", rowNames = F, colNames = T, detectDates = T)

# Importing Gas Chromatpgraphy Mass Spectometry data 
GCMS <-  as.data.frame(read.xlsx("Data/Profiles.xlsx", sheet = "GCMS", rowNames = T, colNames = T))

# Importing Isotopic Ratio Mass Spectrometry data
IRMS <-  as.data.frame(read.xlsx("Data/Profiles.xlsx", sheet = "IRMS", rowNames = T, colNames = T))

# Importing Capillary Electrophoresis data
CE <- as.data.frame(read.xlsx("Data/Profiles.xlsx", sheet = "CE", rowNames = T, colNames = T))



# Choosing target variables - GCMS profiles only: -------------------------

# (1) Presence of variables in all specimens
Count <- data.frame(Variables = colnames(GCMS),
                    Percentage = round((colSums(GCMS !=0)/nrow(GCMS)), digits = 5),
                    row.names = NULL)

Bar <- ggplot(Count, aes(Variables, Percentage*100)) +
  geom_bar(stat = "identity") +
  theme_light()+
  theme(axis.title.x = element_text(face="bold",size = 12,margin=margin(10,0,0,0)),
        axis.title.y = element_text(face="bold",size = 12,margin=margin(0,10,0,0)),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
        panel.grid = element_blank(),
        panel.background=element_rect(fill=NA),
        legend.position = "none")+
  scale_y_continuous("Percentage of specimens (%)")  +
  coord_cartesian(ylim = c(0,100))
Bar

# (2) Intra-varaibility and inter-variability of variables
Specimen_List <- data.frame(row.names = rownames(GCMS))
Specimen_List$Group <- Lookup$Group[match(row.names(Specimen_List),Lookup$Specimen)]
Specimen_List$Count <- plyr::count(Specimen_List,"Group")$freq[match(Specimen_List$Group,plyr::count(Specimen_List,"Group")$Group)]

Lkd <- lapply(c("G_158","G_030","G_162","G_285"),function(SG){
  X = stack(subset(GCMS, row.names(GCMS) %in% row.names(Specimen_List[Specimen_List$Group==SG,])))
  X$Group = rep(SG,nrow(X))
  return(X)
})
ULkd <-  stack(subset(GCMS, row.names(GCMS) %in% row.names(Specimen_List[Specimen_List$Count==1,])))
ULkd$Group <- rep("Multiple",nrow(ULkd))
Group.labs <- c("Inter-variability",
                "Intra-variability of Group A",
                "Intra-variability of Group B",
                "Intra-variability of Group C",
                "Intra-variability of Group D")
names(Group.labs) <- c("Multiple","G_158","G_030","G_162","G_285")

BoxP <- rbind(do.call(rbind,Lkd),ULkd)

psplit <- ggplot(BoxP, aes(x=ind,y=values)) + 
  geom_boxplot(outlier.shape=1,outlier.size=1,lwd=0.3) +
  facet_wrap(~Group, scale= "fixed", ncol=1, labeller = labeller(Group=Group.labs)) +
  labs(x="Target variables",y="Pre-treated variable intensity") +
  theme_light()+
  theme(axis.title.x = element_text(face="bold",size = 12,margin=margin(10,0,0,0)),
        axis.title.y = element_text(face="bold",size = 12,margin=margin(0,10,0,0)),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 0, vjust = 0.5),
        strip.text=element_text(face="bold"),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
psplit

# (3) Correlation of variables
M = cor(GCMS,method = "spearman")
corrplot <- ggcorrplot(M, type = "upper",outline.col = "white", p.mat = ggcorrplot::cor_pmat(M),
                       colors = c("#595959", "white", "#595959"))+
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
        panel.grid.minor = element_blank())
corrplot



# Optimisation of the linked and unlinked specimen populations: -----------

# Reading in several functions needed for subsequent analysis
 source("Files/DataOptFunctions.R", local = T)

#### GCMS 
# Applying population rules and comparison metrics to the data
GCMS_CM_R = 
  lapply(
    list(CAN=CAN,EUC=EUC,MAN=MAN,MCF=MCF,PCC=PCC),
    function(CMs){
      list(
        R1 = R1(CMs(as.matrix(GCMS))),
        R2 = data.frame(R2(CMs(L(GCMS)),CMs(U(GCMS)))),
        R3 = R3(CMs(as.matrix(GCMS))),
        R4 = data.frame(R4(CMs(L(GCMS)),CMs(as.matrix(GCMS)))))
    })

# Adding in comparison metric and pupulation rule info
for(i in 1:length(GCMS_CM_R)){
  for(j in 1:length(GCMS_CM_R[[1]])){
    GCMS_CM_R[[i]][[j]]$Freq = (GCMS_CM_R[[i]][[j]]$Freq/max(GCMS_CM_R[[i]][[j]]$Freq))*100
    GCMS_CM_R[[i]][[j]]$CM = rep(c("CAN","EUC","MAN","MCF","PCC")[i],nrow(GCMS_CM_R[[i]][[j]]))
    GCMS_CM_R[[i]][[j]]$R = rep((paste0("R",j)),nrow(GCMS_CM_R[[i]][[j]])) 
  }
}

# Extracting the area under curve (AUC) of each reciever operating characteristic (ROC) calculation
GCMS_CM_R <- do.call(rbind,do.call(rbind,GCMS_CM_R))

GCMS_AUC <- GCMS_CM_R %>% 
  group_by(CM,R) %>% 
  group_map(~{roc(.x$label, .x$Freq, levels = c("Inter", "Intra"), direction = ">")$auc})

GCMS_ROC <- data.frame(AUC = unlist(GCMS_AUC)) %>% bind_cols(GCMS_CM_R %>% distinct(CM,R))

# Defining the optimal comparison metric and population rule based on the ROC AUC
OPT_GCMS_AUC <- max(GCMS_ROC$AUC)
OPT_GCMS_CM <- GCMS_ROC[which.max(GCMS_ROC$AUC),"CM"]
OPT_GCMS_PR <- GCMS_ROC[which.max(GCMS_ROC$AUC),"R"]

OPT_GCMS_CM_R <- GCMS_CM_R %>%
  filter(CM == OPT_GCMS_CM,
         R == OPT_GCMS_PR)

# Visualising the seperation between linked and unlinked populations for each CM and PR combination
ggplot(GCMS_CM_R, aes(x=label, y=Freq)) + 
  geom_boxplot(outlier.shape=1,outlier.size=1) + 
  facet_grid(CM~R, scales="free")+
  labs(y="CM score")+
  scale_fill_manual(values=c("white","white"))+
  theme_light()+
  theme(legend.position="none",
        axis.title.y=element_text(face="bold"), 
        axis.title.x=element_blank(),
        strip.text=element_text(face="bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

# Density plot of the optimal CM and PR combination
ggplot(OPT_GCMS_CM_R,aes(x=Freq,colour=label))+
  geom_density()+
  coord_cartesian(xlim = c(0,100))+
  labs(x="CM score",y="Frequency (%)")+                          
  scale_color_grey()+
  theme_light()+
  theme(plot.title=element_text(face="bold",hjust=.5),
        axis.title=element_text(face="bold"),
        legend.title=element_blank(),
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())



#### IRMS
# Adjusting neagtive values and applying pre-treatments
IRMS = IRMS-(min(IRMS))+0.1

IRMS_N = IRMS/rowSums(IRMS)
IRMS_N2R = IRMS_N^0.5
IRMS_N4R = IRMS_N^0.25
IRMS_L = log10(IRMS)

# Applying population rules and comparison metrics to the data
IRMS_PT_CM_R = lapply(  
  list(N=IRMS_N,N2R=IRMS_N2R,N4R=IRMS_N4R,L=IRMS_L), 
  function(X){
    lapply(
      list(CAN=CAN,EUC=EUC,MAN=MAN,MCF=MCF,PCC=PCC),
      function(CM,PT){
        list(
          R1 = R1(CM(as.matrix(PT))),
          R2 = data.frame(R2(CM(L(PT)),CM(U(PT)))),
          R3 = R3(CM(as.matrix(PT))),
          R4 = data.frame(R4(CM(L(PT)),CM(as.matrix(PT)))))
      },
      PT=X)
  })

# Adding in comparison metric and pupulation rule info
for(i in 1:length(IRMS_PT_CM_R)){
  for(j in 1:length(IRMS_PT_CM_R[[1]])){
    for(k in 1:length(IRMS_PT_CM_R[[1]][[1]])){
      IRMS_PT_CM_R[[i]][[j]][[k]]$Freq = (IRMS_PT_CM_R[[i]][[j]][[k]]$Freq/max(IRMS_PT_CM_R[[i]][[j]][[k]]$Freq))*100
      IRMS_PT_CM_R[[i]][[j]][[k]]$PT = rep(c("N","N2R","N4R","L")[i],nrow(IRMS_PT_CM_R[[i]][[j]][[k]]))
      IRMS_PT_CM_R[[i]][[j]][[k]]$CM = rep(c("CAN","EUC","MAN","MCF","PCC")[j],nrow(IRMS_PT_CM_R[[i]][[j]][[k]]))
      IRMS_PT_CM_R[[i]][[j]][[k]]$R = rep((paste0("R",k)),nrow(IRMS_PT_CM_R[[i]][[j]][[k]])) 
    }
  }
}

# Extracting the area under curve (AUC) of each reciever operating characteristic (ROC) calculation
IRMS_PT_CM_R <- do.call(rbind,do.call(rbind,do.call(rbind,IRMS_PT_CM_R)))

IRMS_AUC <- IRMS_PT_CM_R %>% 
  group_by(PT,CM,R) %>% 
  group_map(~{roc(.x$label, .x$Freq, levels = c("Inter", "Intra"), direction = ">")$auc})

IRMS_ROC <- data.frame(AUC = unlist(IRMS_AUC)) %>% bind_cols(IRMS_PT_CM_R %>% distinct(PT,CM,R))

# Defining the optimal comparison metric and population rule based on the ROC AUC
OPT_IRMS_AUC <- max(IRMS_ROC$AUC)
OPT_IRMS_PT <- IRMS_ROC[which.max(IRMS_ROC$AUC),"PT"]
OPT_IRMS_CM <- IRMS_ROC[which.max(IRMS_ROC$AUC),"CM"]
OPT_IRMS_PR <- IRMS_ROC[which.max(IRMS_ROC$AUC),"R"]

OPT_IRMS_PT_CM_R <- IRMS_PT_CM_R %>%
  filter(PT == OPT_IRMS_PT,
         CM == OPT_IRMS_CM,
         R == OPT_IRMS_PR)

# Visualising the seperation between linked and unlinked populations for each CM and PR combination
ggplot(IRMS_PT_CM_R, aes(x=label, y=Freq)) + 
  geom_boxplot(outlier.shape=1,outlier.size=1) + 
  facet_nested(PT+CM~R, scales="free")+
  labs(y="CM score")+
  scale_fill_manual(values=c("white","white"))+
  theme_light()+
  theme(legend.position="none",
        axis.title.y=element_text(face="bold"), 
        axis.title.x=element_blank(),
        strip.text=element_text(face="bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

# Density plot of the optimal CM and PR combination
ggplot(OPT_IRMS_PT_CM_R,aes(x=Freq,colour=label))+
  geom_density()+
  coord_cartesian(xlim = c(0,100))+
  labs(x="CM score",y="Frequency (%)")+
  scale_color_grey()+
  theme_light()+
  theme(plot.title=element_text(face="bold",hjust=.5),
        axis.title=element_text(face="bold"),
        legend.title=element_blank(),
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())



#### CE
# Applying pre-treatments
CE_N = CE/rowSums(CE)
CE_N2R = CE_N^0.5
CE_N4R = CE_N^0.25

# Applying population rules and comparison metrics to the data
CE_PT_CM_R = lapply(  
  list(N=CE_N,N2R=CE_N2R,N4R=CE_N4R), 
  function(X){
    lapply(
      list(CAN=CAN,EUC=EUC,MAN=MAN,MCF=MCF,PCC=PCC),
      function(CM,PT){
        list(
          R1 = R1(CM(as.matrix(PT))),
          R2 = data.frame(R2(CM(L(PT)),CM(U(PT)))),
          R3 = R3(CM(as.matrix(PT))),
          R4 = data.frame(R4(CM(L(PT)),CM(as.matrix(PT)))))
      },
      PT=X)
  })

# Adding in comparison metric and pupulation rule info
for(i in 1:length(CE_PT_CM_R)){
  for(j in 1:length(CE_PT_CM_R[[1]])){
    for(k in 1:length(CE_PT_CM_R[[1]][[1]])){
      CE_PT_CM_R[[i]][[j]][[k]]$Freq = (CE_PT_CM_R[[i]][[j]][[k]]$Freq/max(CE_PT_CM_R[[i]][[j]][[k]]$Freq))*100
      CE_PT_CM_R[[i]][[j]][[k]]$PT = rep(c("N","N2R","N4R")[i],nrow(CE_PT_CM_R[[i]][[j]][[k]]))
      CE_PT_CM_R[[i]][[j]][[k]]$CM = rep(c("CAN","EUC","MAN","MCF","PCC")[j],nrow(CE_PT_CM_R[[i]][[j]][[k]]))
      CE_PT_CM_R[[i]][[j]][[k]]$R = rep((paste0("R",k)),nrow(CE_PT_CM_R[[i]][[j]][[k]])) 
    }
  }
}

# Extracting the area under curve (AUC) of each reciever operating characteristic (ROC) calculation
CE_PT_CM_R <- do.call(rbind,do.call(rbind,do.call(rbind,CE_PT_CM_R)))

CE_AUC <- CE_PT_CM_R %>% 
  group_by(PT,CM,R) %>% 
  group_map(~{roc(.x$label, .x$Freq, levels = c("Inter", "Intra"), direction = ">")$auc})

CE_ROC <- data.frame(AUC = unlist(CE_AUC)) %>% bind_cols(CE_PT_CM_R %>% distinct(PT,CM,R))

# Defining the optimal comparison metric and population rule based on the ROC AUC
OPT_CE_AUC <- max(CE_ROC$AUC)
OPT_CE_PT <- CE_ROC[which.max(CE_ROC$AUC),"PT"]
OPT_CE_CM <- CE_ROC[which.max(CE_ROC$AUC),"CM"]
OPT_CE_PR <- CE_ROC[which.max(CE_ROC$AUC),"R"]

OPT_CE_PT_CM_R <- CE_PT_CM_R %>%
  filter(PT == OPT_CE_PT,
         CM == OPT_CE_CM,
         R == OPT_CE_PR)

# Visualising the seperation between linked and unlinked populations for each CM and PR combination
ggplot(CE_PT_CM_R, aes(x=label, y=Freq)) + 
  geom_boxplot(outlier.shape=1,outlier.size=1) + 
  facet_nested(PT+CM~R, scales="free")+
  labs(y="CM score")+
  scale_fill_manual(values=c("white","white"))+
  theme_light()+
  theme(legend.position="none",
        axis.title.y=element_text(face="bold"), 
        axis.title.x=element_blank(),
        strip.text=element_text(face="bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

# Density plot of the optimal CM and PR combination
ggplot(OPT_CE_PT_CM_R,aes(x=Freq,colour=label))+
  geom_density()+
  coord_cartesian(xlim = c(0,100))+
  labs(x="CM score",y="Frequency (%)")+
  scale_color_grey()+
  theme_light()+
  theme(plot.title=element_text(face="bold",hjust=.5),
        axis.title=element_text(face="bold"),
        legend.title=element_blank(),
        legend.position="none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())




# Defining the optimal analytical technique (i.e. GCMS, IRMS or CE) -------
OPT_AT <- which.max(c(GCMS = OPT_GCMS_AUC, IRMS = OPT_IRMS_AUC, CE = OPT_CE_AUC))

# Concluding statement
paste0("The dataset with the largest discrimination between the linked and unlinked populations is the ",names(OPT_AT),
       " dataset. In this instance the optimal combination of comparison metric and population rule is ",
       get(paste0("OPT_",names(OPT_AT),"_CM")), " and ", get(paste0("OPT_",names(OPT_AT),"_PR")), ", respectively.")