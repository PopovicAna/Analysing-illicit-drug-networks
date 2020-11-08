#Linked (L) and Unlinked (U) population definitions
L <- function(x){as.matrix(subset(x, row.names(x) %in% row.names(Specimen_List[Specimen_List$Count>1 & Specimen_List$Count<11,])))}
U <- function(x){as.matrix(subset(x, row.names(x) %in% row.names(Specimen_List[Specimen_List$Count==1,])))}

#Population rules (PR) 
R1 = function(CM){
  CM[upper.tri(CM)] = NA
  diag(CM) = NA
  CM = as.data.frame.table(CM)
  CM = na.omit(CM)
  CM$Group1 = Lookup$Group[match(CM$Var1, Lookup$Specimen)]
  CM$Group2 = Lookup$Group[match(CM$Var2, Lookup$Specimen)]
  CM$label = ifelse(CM$Group1 == CM$Group2, "Intra", "Inter")
  CM=CM[,c("Freq","label")]
}
R2 = function(CML,CMU){
  CML[upper.tri(CML)] = NA
  diag(CML) = NA
  CML = as.data.frame.table(CML)
  CML = na.omit(CML)
  CML$Group1 = Lookup$Group[match(CML$Var1, Lookup$Specimen)]
  CML$Group2 = Lookup$Group[match(CML$Var2, Lookup$Specimen)]
  CML$label = ifelse(CML$Group1 == CML$Group2, "Intra", "Inter")
  CML = subset(CML, !(CML$label == "Inter"))
  CMU[upper.tri(CMU)] = NA
  diag(CMU) = NA
  CMU = as.data.frame.table(CMU)
  CMU = na.omit(CMU)
  CMU$Group1 = Lookup$Group[match(CMU$Var1, Lookup$Specimen)]
  CMU$Group2 = Lookup$Group[match(CMU$Var2, Lookup$Specimen)]
  CMU$label = ifelse(CMU$Group1 == CMU$Group2, "Intra", "Inter")
  CMU = subset(CMU, !(CMU$label == "Intra"))
  rbind(CML[,c("Freq","label")],CMU[,c("Freq","label")])
}
R3 = function(CM){
  CM[upper.tri(CM)] = NA
  diag(CM) = NA
  CM = as.data.frame.table(CM)
  CM = na.omit(CM)
  CM$Group1 = Lookup$Group[match(CM$Var1, Lookup$Specimen)]
  CM$Group2 = Lookup$Group[match(CM$Var2, Lookup$Specimen)]
  CM$PriR1 = Lookup$PriR[match(CM$Var1, Lookup$Specimen)]
  CM$PriR2 = Lookup$PriR[match(CM$Var2, Lookup$Specimen)]
  CM$Pre1 = Lookup$Pre[match(CM$Var1, Lookup$Specimen)]
  CM$Pre2 = Lookup$Pre[match(CM$Var2, Lookup$Specimen)]
  CM$PreR1 = Lookup$PreR[match(CM$Var1, Lookup$Specimen)]
  CM$PreR2 = Lookup$PreR[match(CM$Var2, Lookup$Specimen)]
  CM$label = ifelse(CM$Group1 == CM$Group2 & CM$PriR1 == CM$PriR2 
                    & CM$Pre1 == CM$Pre2 & CM$PreR1 == CM$PreR2, "Intra",
                    ifelse(CM$Group1 == CM$Group2 | CM$PriR1 == CM$PriR2 
                           | CM$Pre1 == CM$Pre2 | CM$PreR1 == CM$PreR2, NA, "Inter"))
  CM = na.omit(CM)
  CM = CM[,c("Freq","label")]
}

R4 = function(CML,CMU){
  CML[upper.tri(CML)] = NA
  diag(CML) = NA
  CML = as.data.frame.table(CML)
  CML = na.omit(CML)
  CML$Group1 = Lookup$Group[match(CML$Var1, Lookup$Specimen)]
  CML$Group2 = Lookup$Group[match(CML$Var2, Lookup$Specimen)]
  CML$PriR1 = Lookup$PriR[match(CML$Var1, Lookup$Specimen)]
  CML$PriR2 = Lookup$PriR[match(CML$Var2, Lookup$Specimen)]
  CML$Pre1 = Lookup$Pre[match(CML$Var1, Lookup$Specimen)]
  CML$Pre2 = Lookup$Pre[match(CML$Var2, Lookup$Specimen)]
  CML$PreR1 = Lookup$PreR[match(CML$Var1, Lookup$Specimen)]
  CML$PreR2 = Lookup$PreR[match(CML$Var2, Lookup$Specimen)]
  CML$label = ifelse(CML$Group1 == CML$Group2 & CML$PriR1 == CML$PriR2 
                     & CML$Pre1 == CML$Pre2 & CML$PreR1 == CML$PreR2, "Intra", "Inter")
  CML = subset(CML, !(CML$label == "Inter"))
  
  CMU[upper.tri(CMU)] = NA
  diag(CMU) = NA
  CMU = as.data.frame.table(CMU)
  CMU = na.omit(CMU)
  CMU$Group1 = Lookup$Group[match(CMU$Var1, Lookup$Specimen)]
  CMU$Group2 = Lookup$Group[match(CMU$Var2, Lookup$Specimen)]
  CMU$PriR1 = Lookup$PriR[match(CMU$Var1, Lookup$Specimen)]
  CMU$PriR2 = Lookup$PriR[match(CMU$Var2, Lookup$Specimen)]
  CMU$Pre1 = Lookup$Pre[match(CMU$Var1, Lookup$Specimen)]
  CMU$Pre2 = Lookup$Pre[match(CMU$Var2, Lookup$Specimen)]
  CMU$PreR1 = Lookup$PreR[match(CMU$Var1, Lookup$Specimen)]
  CMU$PreR2 = Lookup$PreR[match(CMU$Var2, Lookup$Specimen)]
  CMU$label = ifelse(CMU$Group1 == CMU$Group2 | CMU$PriR1 == CMU$PriR2 |
                       CMU$Pre1 == CMU$Pre2 | CMU$PreR1 == CMU$PreR2,"Intra", "Inter")
  CMU = subset(CMU, !(CMU$label == "Intra"))
  
  rbind(CML[,c("Freq","label")],CMU[,c("Freq","label")])
}

#Comparison metrics
MCF <- function(V){
  (100-(100*((((V)%*%(t(V)))^2)/(rowSums((V)^2))%*%(t(rowSums((V)^2))))))
}
PCC = function(V){
  (((1-(cor(t(V), method = "pearson")))/2)*100)
}
EUC = function(V){
  as.matrix(dist(V, method = "euclidean"))
}
MAN = function(V){
  as.matrix(dist(V, method = "manhattan"))
}
CAN = function(V){
  as.matrix(dist(V, method = "canberra"))
  }

#ROC AUC
RCM = function(CM){
  roc(CM$label, CM$Freq, levels = c("Inter", "Intra"), direction = ">")
}