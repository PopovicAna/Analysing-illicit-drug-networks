library(igraph)

source("Files/Prioritisation-of-analytical-techniques.R", local = T)


# DEEP LEVEL --------------------------------------------------------------

# Deterministic approach: -------------------------------------------------

# Defining TP, TN, FP and FN rates at each threshold values (THV - aka CM score)
Cutoffs = data.frame(
  THV = RCM(OPT_GCMS_CM_R)$thresholds, 
  TPR = RCM(OPT_GCMS_CM_R)$sensitivities, 
  FPR = 1-RCM(OPT_GCMS_CM_R)$specificities,
  FNR = 1-RCM(OPT_GCMS_CM_R)$sensitivities, 
  TNR = RCM(OPT_GCMS_CM_R)$specificities)

# Visualising the FPR for each THV
ggplot(Cutoffs)+
  geom_line(aes(x=THV,y=FPR),size=1,lty=1)+
  geom_line(aes(x=THV,y=FNR),size=1,lty=2)+
  scale_fill_grey()+
  theme_light()+
  labs(x="CM score",y="Ratio")+
  theme(axis.title=element_text(face="bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())      

# Setting an acceptable FPR to define whether specimen pairs are linked or not
# Note: For the purpose of this research an acceptable FPR is 0.025 (i.e. 2.5%)
LINK_THV <- as.numeric(tail(subset(Cutoffs, FPR > 0.025, select = THV), n = 1))



# Bayesian approach: ------------------------------------------------------

# Density plot of the optimal CM and PR combination 
OPT_DENS <- ggplot(OPT_GCMS_CM_R,aes(x=Freq,colour=label))+
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
OPT_DENS

# Extracting the frequency (y-axis) of the linked and unlinked populations at each THV (x-axis)
GG_OPT_DENS = ggplot_build(OPT_DENS)

# Defining the likelihood ratio (LR) of a link at each THV
LRs = data.frame(THV=GG_OPT_DENS$data[[1]]$x[1:512], 
                 LR=((GG_OPT_DENS$data[[1]]$y[GG_OPT_DENS$data[[1]]$group==2])/
                       (GG_OPT_DENS$data[[1]]$y[GG_OPT_DENS$data[[1]]$group==1])))

# Visualising the LR for each THV
ggplot(LRs, aes(x=THV,y=LR))+
  geom_line(size=1)+
  scale_fill_grey()+
  theme_light()+
  labs(x="CM score",y="LR")+
  theme(axis.title=element_text(face="bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

# Function for the verbal equivalent of the LR
LR_Verbal <- function(LR){
  ifelse(LR>10000,"VSS_H1",
         ifelse(LR<=10000&LR>1000,"SS_H1",
                ifelse(LR<=1000&LR>100,"MSS_H1",
                       ifelse(LR<=100&LR>10,"MS_H1",
                              ifelse(LR<=10&LR>1,"LS_H1",
                                     ifelse(LR<=1&LR>0.1,"LS_H2",
                                            ifelse(LR<=0.1&LR>0.01,"MS_H2",
                                                   ifelse(LR<=0.01&LR>0.001,"MSS_H2",
                                                          ifelse(LR<=0.001&LR>0.0001,"SS_H2",
                                                                 ifelse(LR<=0.0001,"VSS_H2","")
                                                          )))))))))
}




# WORKING LEVEL -----------------------------------------------------------

# Extracting a subset of values
Sub_GCMS <- GCMS %>% 
  filter(rownames(GCMS) %in% Lookup$Specimen[Lookup$Date=="2017-09-01"])
rownames(Sub_GCMS) = Lookup$Specimen[Lookup$Date=="2017-09-01"]

# Calculating similarity between specimens based on the optimal comparison metric
Scores_Sub_GCMS <- get(OPT_GCMS_CM)(Sub_GCMS)
Scores_Sub_GCMS[upper.tri(Scores_Sub_GCMS)] = NA
diag(Scores_Sub_GCMS) = NA
Scores_Sub_GCMS = as.data.frame.table(Scores_Sub_GCMS)
Scores_Sub_GCMS = na.omit(Scores_Sub_GCMS)
colnames(Scores_Sub_GCMS) <- c("From","To","Score")

#Defining the LR and verbal equivalent based on the deep level results
Scores_Sub_GCMS$LR <-  unlist(lapply(Scores_Sub_GCMS$Score,function(W) head(subset(LRs, THV>= W, select = LR),n=1)))
Scores_Sub_GCMS$LR_V <- unlist(LR_Verbal(Scores_Sub_GCMS$LR))

# Plot of links based on the Deterministic approach
netTHV = graph_from_data_frame(d=Scores_Sub_GCMS,directed = F)
E(netTHV)$lty="solid"
E(netTHV)$width="2"
E(netTHV)$color="darkgrey"

LayOut <- matrix(c(2,0,2,0,0,2,1,2,1,0,
                   3,2.5,2.5,0.5,1.5,1.5,1,0.5,0,3),
                 nrow = 10,
                 ncol = 2)

# plot(delete_edges(netTHV, which(E(netTHV)$Score>=LINK_THV)),layout = LayOut)

# Plot of the likelihood of a link based on the Bayesian approach
set.seed(1272)
netLR = graph_from_data_frame(d=Scores_Sub_GCMS,directed = F)
E(netLR)$lty=
  ifelse(Scores_Sub_GCMS$LR>1000,"solid",
         ifelse(Scores_Sub_GCMS$LR<=1000&Scores_Sub_GCMS$LR>100,"dashed",
                ifelse(Scores_Sub_GCMS$LR<=100&Scores_Sub_GCMS$LR>10,"dotdash",
                       ifelse(Scores_Sub_GCMS$LR<=10&Scores_Sub_GCMS$LR>1,"dotted",
                              ifelse(Scores_Sub_GCMS$LR<=1&Scores_Sub_GCMS$LR>0.1,"dotted",
                                     ifelse(Scores_Sub_GCMS$LR<=0.1&Scores_Sub_GCMS$LR>0.01,"dotdash",
                                            ifelse(Scores_Sub_GCMS$LR<=0.01&Scores_Sub_GCMS$LR>0.001,"dashed",
                                                   ifelse(Scores_Sub_GCMS$LR<=0.001&Scores_Sub_GCMS$LR>0.0001,"solid",
                                                          ifelse(Scores_Sub_GCMS$LR<=0.0001,"solid","")
                                                   ))))))))
E(netLR)$width=ifelse(Scores_Sub_GCMS$LR<0.0001,"4","2")
E(netLR)$color=ifelse(Scores_Sub_GCMS$LR>1,"Black",ifelse(Scores_Sub_GCMS$LR<1,"Red",""))
# plot(netLR,layout = LayOut)

# Remove any edges with no support for linkage between specimens
#plot(delete_edges(netLR, which(E(netLR)$LR<=1)),layout = LayOut))

# Visualising the Det and Bayes results for the subset of specimens -------
par(mfrow=c(1,2))
plot(delete_edges(netTHV, which(E(netTHV)$Score>=LINK_THV)), layout = LayOut)
plot(delete_edges(netLR, which(E(netLR)$LR<=1)), layout = LayOut)




"Make hypothesis based on the connections seen"