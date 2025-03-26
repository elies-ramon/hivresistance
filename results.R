####################
##### RESULTS ######
####################

# Llibreries #
library(kerntools)
library(cowplot)
library(ggplot2)
library(reshape2)


# Declaració de funcions i constants #
source("helper_functions.R")


# Càrrega dels indexes #
load("Partitions/PI_IDX")
N <- length(c(IDX$work,IDX$additional_test_subtype))
N

# Càrrega de les posicions de resistència #

WTSeqseq <-toSeq(WTSeq)

# Càrrega de les dades #
load("Data/PI_MultiA_clean")
rownames(DATA) <- 1:nrow(DATA)

te_add_subtype <-  DATA[IDX$additional_test_subtype,"SeqID"]
te_idx <- DATA[ IDX$test,"SeqID"]
anyDuplicated(te_idx)
anyDuplicated(te_add_subtype)

DATA <- droplevels(DATA[as.numeric(c(IDX$work,IDX$additional_test_subtype)),])
dim(DATA)  #2212  115

DATAseq <- toSeq(DATA[,PrPositions])

Partition <- rep("Training",nrow(DATA))
Partition[which( DATA$SeqID %in% te_idx)] <- "Test"
Partition[which( DATA$SeqID %in% te_add_subtype)] <- "Additional test"
summary(factor(Partition))
sum(DATA[Partition=="Additional test","Subtype"] == "B")

load(file="Kernel_matrices/PI_Kint.dat")
load(file="Kernel_matrices/PI_mixSP1.dat")
load(file="Kernel_matrices/PI_mixSP2.dat")
load(file="Kernel_matrices/PI_mixSP3.dat")
load(file="Kernel_matrices/PI_MIXLONG1.dat")



# Figures #

### Histogram of kernel matrices

par(mfrow=c(3,2))
histK(cosNorm(mixSP1$K),vn = TRUE,main = "Histogram of Spectrum l=1",freq = FALSE ,border="white",col="grey60")
histK(cosNorm(Kint$K),vn = TRUE,main = "Histogram of intersect",freq = FALSE ,border="white",col="grey60")

histK(cosNorm(mixSP2$K),vn = TRUE,main = "Histogram of Spectrum l=2",freq = FALSE ,border="white",col="grey60")
histK(cosNorm(Mixlong1$K),vn = TRUE,main = "Histogram of interactions",freq = FALSE ,border="white",col="grey60")

histK(cosNorm(mixSP3$K),vn = TRUE,main = "Histogram of Spectrum l=3",freq = FALSE ,border="white",col="grey60")
par(mfrow=c(1,1))


#### KPCAs
drug <- "LPV"

load(file="Results/PI_PCS_int.dat")
load(file="Results/PI_PCS_mixSP1.dat")
load(file="Results/PI_PCS_mixSP2.dat")
load(file="Results/PI_PCS_mixSP3.dat")

allele_names <- paste0(sub("^P","",sub('_.*', '',  colnames(DATAdum_delete))),sub('.*_', '', colnames(DATAdum_delete)))
colnames(pcs_int$loadings) <- allele_names

pca_int <- kPCA(Kint$K,center=TRUE,plot=1:2,colors = "grey99")
pca_sp1 <- kPCA(mixSP1$K,center=TRUE,plot=1:2,colors = "grey99")
pca_sp2 <- kPCA(mixSP2$K,center=TRUE,plot=1:2,colors = "grey99")
pca_sp3 <- kPCA(mixSP3$K,center=TRUE,plot=1:2,colors = "grey99")
pca_long <- kPCA(Mixlong1$K,center=TRUE,plot=1:2,colors = "grey99")

### Figura SP1

partitioncolors <- c("#6DB6FF","#996035","#FFB6DB")
panelA <-   pca_sp1$plot + ggplot2::ggtitle("Spectrum k=1: Partitions") +
  ggplot2::geom_point(ggplot2::aes(colour =  Partition)) + scale_colour_manual(values=partitioncolors) + theme(legend.position = "bottom")

panelB <- pca_sp1$plot + ggplot2::ggtitle(paste("Spectrum k=1:",drug)) +
  ggplot2::geom_point(ggplot2::aes(colour = DATA[,drug])) + theme(legend.position = "bottom") +
  ggplot2::scale_colour_gradientn(name=drug,colours=c("green","yellow","red"), na.value = "grey70")
panelB <- kPCA_arrows(panelB, contributions = t(pcs_s1$loadings[1:2,c("V","I","L","R","K","G")]),
            colour = "grey30")

plot_grid(panelB,panelA,ncol =2,labels=c("A","B"))


### Figura Partitions

panelA <-   pca_int$plot + ggplot2::ggtitle("Intersect kPCA - Partitions") +
  ggplot2::geom_point(ggplot2::aes(colour =  Partition),show.legend = FALSE) + scale_colour_manual(values=partitioncolors)

panelB <-   pca_sp2$plot + ggplot2::ggtitle("Spectrum k=2 kPCA - Partitions") +
  ggplot2::geom_point(ggplot2::aes(colour =  Partition),show.legend = FALSE) + scale_colour_manual(values=partitioncolors)
panelC <-   pca_long$plot + ggplot2::ggtitle("Interactions kPCA - Partitions") +
  ggplot2::geom_point(ggplot2::aes(colour =  Partition)) + scale_colour_manual(values=partitioncolors)
panelD <- pca_sp3$plot + ggplot2::ggtitle("Spectrum k=3 kPCA - Partitions") +
  ggplot2::geom_point(ggplot2::aes(colour =  Partition),show.legend = FALSE) +scale_colour_manual(values=partitioncolors)

legend <- ggplotGrob(panelC + theme(legend.position = "bottom"))$grobs
legend <- legend[[which(sapply(legend, function(x) x$name) == "guide-box")]]

 
plot_grid(plot_grid(panelA,panelB,panelC+ ggplot2::theme(legend.position="none"),panelD,nrow=2,ncol=2,
          labels = LETTERS[1:4]), legend,nrow=2,rel_heights = c(0.9,0.1))


#sutipus del HIV-1
subtype <- DATA[,"Subtype"]
subtype[subtype=="B"] <- NA
kPCA(Kint$K,center=TRUE,plot=1:2,y =subtype)
kPCA(mixSP3$K,center=TRUE,plot=1:2,y =subtype)
 

### Figura PCs

color1 <- color2 <- color3 <- color4 <- rep("#F2DACD",20)
color1[which(names(sort(abs(pcs_int$loadings[1,]),decreasing = T)[1:20]) %in% paste0(1:99,WTSeq))] <-  "lightskyblue1"
color2[vect_grep(names(sort(abs(pcs_s2$loadings[1,]),decreasing = T)[1:20]),WTSeqseq)] <-  "lightskyblue1"
color3[vect_grep(names(sort(abs(pcs_s3$loadings[1,]),decreasing = T)[1:20]),WTSeqseq)] <-  "lightskyblue1"
color4[vect_grep(names(sort(abs(pcs_s3$loadings[2,]),decreasing = T)[1:20]),WTSeqseq)] <-  "lightskyblue1"


par(mfrow=c(2,2))
pc1_list <- plotImp(pcs_int$loadings[1,], absolute=FALSE, relative = FALSE,  main="Intersect PC1",
                    leftmargin = 5,  color=rev(color1), nfeat=20)

pc1_list2 <- plotImp(pcs_s2$loadings[1,], absolute=FALSE, relative = FALSE,  main="Spectrum k=2 PC1",
                     leftmargin = 5,  color=rev(color2), nfeat=20)

pc1_list3 <- plotImp(pcs_s3$loadings[1,], absolute=FALSE, relative = FALSE,  main="Spectrum k=3 PC1",
                     leftmargin = 5,  color=rev(color3), nfeat=20)
pc2_list3 <- plotImp(pcs_s3$loadings[2,], absolute=FALSE, relative = FALSE,  main="Spectrum k=3 PC2",
                     leftmargin = 5, color=rev(color4), nfeat=20)

pc2_list <- plotImp(pcs_int$loadings[2,], absolute=FALSE, relative = FALSE,  main="PC2",
                    leftmargin = 5,  color="pink",nfeat=20)

pc2_list2 <- plotImp(pcs_s2$loadings[2,], absolute=FALSE, relative = FALSE,  main="PC2",
                     leftmargin = 5,  color="pink",nfeat=20)

### Figura LPV

panelA <- pca_int$plot + ggplot2::ggtitle(paste("Intersect kPCA -",drug)) +
  ggplot2::geom_point(ggplot2::aes(colour = DATA[,drug]),show.legend = FALSE) +
  ggplot2::scale_colour_gradientn(name=drug,colours=c("green","yellow","red"), na.value = "grey70")
panelA <- kPCA_arrows(panelA,
                  contributions =  t(pcs_int$loadings[1:2,union(pc2_list$first_features,pc1_list$first_features)]),
                  colour = "grey40")
panelB <- pca_sp2$plot + ggplot2::ggtitle(paste("Spectrum k=2 kPCA -",drug)) +
  ggplot2::geom_point(ggplot2::aes(colour = DATA[,drug]),show.legend = FALSE) +
  ggplot2::scale_colour_gradientn(name=drug,colours=c("green","yellow","red"), na.value = "grey70")
panelB <- kPCA_arrows(panelB,
                  contributions = t(pcs_s2$loadings[1:2,union(pc2_list2$first_features,pc1_list2$first_features)]),
                  colour = "grey40")
panelC <- pca_long$plot + ggplot2::ggtitle(paste("Interaction kPCA -",drug)) +
  ggplot2::geom_point(ggplot2::aes(colour = DATA[,drug])) +
  ggplot2::scale_colour_gradientn(name=drug,colours=c("green","yellow","red"), na.value = "grey70")
panelD <- pca_sp3$plot + ggplot2::ggtitle(paste("Spectrum k=3 kPCA -",drug)) +
  ggplot2::geom_point(ggplot2::aes(colour = DATA[,drug]),show.legend = FALSE) +
  ggplot2::scale_colour_gradientn(name=drug,colours=c("green","yellow","red"), na.value = "grey70")
panelD <- kPCA_arrows(panelD,
                  contributions = t(pcs_s3$loadings[1:2,union(pc2_list3$first_features[1:18],pc1_list3$first_features[1:8])]),
                  colour = "grey40")

legend <- ggplotGrob(panelC + theme(legend.position = "bottom"))$grobs
legend <- legend[[which(sapply(legend, function(x) x$name) == "guide-box")]]


plot_grid(plot_grid(panelA,panelB,panelC+ggplot2::theme(legend.position="none"),panelD,nrow=2,
          labels = LETTERS[1:4]), legend,nrow=2,rel_heights = c(0.9,0.1))

for(drug in DRUGS[-4]) {
  panelA <- pca_int$plot + ggplot2::ggtitle(paste("Intersect kPCA -",drug)) +
    ggplot2::geom_point(ggplot2::aes(colour = DATA[,drug]),show.legend = FALSE) +
    ggplot2::scale_colour_gradientn(name=drug,colours=c("green","yellow","red"), na.value = "grey70")
  panelB <- pca_sp2$plot + ggplot2::ggtitle(paste("Spectrum k=2 kPCA -",drug)) +
    ggplot2::geom_point(ggplot2::aes(colour = DATA[,drug]),show.legend = FALSE) +
    ggplot2::scale_colour_gradientn(name=drug,colours=c("green","yellow","red"), na.value = "grey70")
  panelC <- pca_long$plot + ggplot2::ggtitle(paste("Interaction kPCA -",drug)) +
    ggplot2::geom_point(ggplot2::aes(colour = DATA[,drug])) +
    ggplot2::scale_colour_gradientn(name=drug,colours=c("green","yellow","red"), na.value = "grey70")
  panelD <- pca_sp3$plot + ggplot2::ggtitle(paste("Spectrum k=3 kPCA -",drug)) +
    ggplot2::geom_point(ggplot2::aes(colour = DATA[,drug]),show.legend = FALSE) +
    ggplot2::scale_colour_gradientn(name=drug,colours=c("green","yellow","red"), na.value = "grey70")

  legend <- ggplotGrob(panelC + theme(legend.position = "bottom"))$grobs
  legend <- legend[[which(sapply(legend, function(x) x$name) == "guide-box")]]

  plot_grid(plot_grid(panelA,panelB,panelC+ggplot2::theme(legend.position="none"),panelD,nrow=2,
                      labels = LETTERS[1:4]), legend,nrow=2,rel_heights = c(0.9,0.1))
}


# Estudi dels PCS #

pcs_s1$loadings[1:2,]
pcs_s2$loadings[1:2,]
pcs_s3$loadings[1:2,]

pc1_list2$first_features

pc1_list3$first_features


sp3_pc2 <- plotImp(pcs_s3$loadings[2,],  nfeat=20,absolute=FALSE, relative = FALSE,
                   main="PC2", leftmargin = 5,  color="pink",y=pcs_s3$loadings[1,])
sp3_pc2$first_features


#### Comparació amb les posicions rellevants per a la resistència

which(!vect_grep(sp2_pc1$first_features,WTSeqseq)) 

which(!vect_grep(sp3_pc1$first_features,WTSeqseq)) 

which(!vect_grep(sp3_pc2$first_features,WTSeqseq))

#  # #### dreta (20 primeres features)
 kmers2 <- stringr::str_locate_all(WTSeqseq, pc1_list2$first_features)
 names(kmers2) <- pc1_list2$first_features
 kmers2
 

data.frame(Positions=  pc1_list3$first_features,stringr::str_locate(WTSeqseq, pc1_list3$first_features))
kmers3 <- stringr::str_locate_all(WTSeqseq, pc1_list3$first_features)
names(kmers3) <- pc1_list3$first_features
kmers3


#### Intersect kernel
pc1_list$first_features [which(pc1_list$first_features %in% paste(names(WTSeq),WTSeq,sep="_"))] 

pc1_list$first_features[which(!(pc1_list$first_features %in% paste(names(WTSeq),WTSeq,sep="_")))] 


#### SVMs

### Performance (validació)
load("Results/val_idx.dat")
load("Results/validation_results.dat")
nmodels <- 5

perform5x2 <- array(NA,dim=c(nmodels,5,length(DRUGS)),dimnames =
                      list(c("sp1","sp2","sp3","kint","mixlong"),1:5,DRUGS))

for(drug in DRUGS)    {
  for(j in 1:nmodels) {
    for(i in 1:5) {
      ok_idx <- val_idx[[i]][complete.cases(DATA[val_idx[[i]],drug])]
      perform5x2[j,i,drug] <-  nmse(DATA[ok_idx,drug],VALS[ok_idx,drug,j])
    }

  }
  print(drug)
  print(rowMeans (perform5x2[ , ,drug]) )
  print(apply (perform5x2[ , ,drug],1,min) )
  print(apply (perform5x2[ , ,drug],1,max) )
}

## Average of the 5 validation errors:
round(t(apply(aperm(perform5x2,c(3,1,2)), 1:2,mean)),digits=3)


### Performance (test)
load("Results/PI_FIT.dat")
load("Results/PI_PREDS.dat")
load("Results/PI_IMPS.dat")
te_idx <-  which(Partition == "Test")

performance <- matrix(0,nrow=length(DRUGS),ncol=6)
rownames(performance) <- DRUGS
colnames(performance)  <- c("Spectrum1","Spectrum2","Spectrum3","Intersect","Interactions","Nte")

for(drug in DRUGS) {
  print(drug)
  Ytest <- DATA[te_idx,drug]
  cc <- complete.cases(Ytest)
  performance[drug,"Nte"] <- sum(cc)
  for(j in 1:5)    performance[drug,j] <- nmse(Ytest[cc],PREDS[te_idx[cc],drug,j])
}

### amb CIs. Comprovar que la predicció és bona en tot el domini
performance_ci <- array(0,dim=c(length(DRUGS),3,5),
                        dimnames =  list(DRUGS,
                                         c("Mean","2.5%","97.5%"),
                                         c("Spectrum1","Spectrum2","Spectrum3","Intersect","Interactions")))
for(drug in DRUGS) {
  print(drug)
  Ytest <- DATA[te_idx,drug]
  cc <- complete.cases(Ytest)
  for(j in 1:5)    performance_ci[drug,,j] <- Boots_CI(target=Ytest[cc], pred=PREDS[te_idx[cc],drug,j],
                                                    index = "nmse", nboots=3000, confidence = 95)
}

#### Taula
PERF <- performance
for(j in 1:5) PERF[,j] <- apply(round(performance_ci[,,j],digits=3),1,function(x)paste0(x[1]," (",x[2],", ",x[3],")"))
PERF
write.csv(PERF,file = "Results/performance_test.csv")


#### Imatge

rownames(perform5x2) <- c("Spectrum1","Spectrum2","Spectrum3","Intersect","Interactions")
box <-  (melt(perform5x2))
colnames(box) <- c("Method","i","Drug","NMSE")

p <- ggplot(box, aes(x = factor(Drug), y = NMSE  , fill = factor(Method))) +
  geom_point(aes(fill = factor(Method)), size = 3, shape = 21, position = position_jitterdodge()) + theme_minimal()+
  xlab("Drug") + ylab("NMSE")+ggtitle("5x2cv")+ylim(c(min(perform5x2,performance_ci),max(perform5x2,performance_ci)))+
  theme(  legend.position = "none") + scale_fill_manual(name="Kernel",values=c("firebrick1","sienna1","goldenrod1","plum","skyblue"))

perf_plot <- data.frame(melt(performance_ci[,1,]),low=melt(performance_ci[,2,])[,3],up=melt(performance_ci[,3,])[,3])

q <- ggplot(perf_plot, aes(x=Var1, y=value, group=Var2, color=Var2)) +
  geom_point(size=1,position = position_jitterdodge())+
  geom_errorbar(aes(ymin=low, ymax=up), width=.5,
                position=position_dodge(0.7)) + theme_minimal()+
  xlab("Drug") + ylab("NMSE")+ggtitle("Test")+ylim(c(min(perform5x2,performance_ci),max(perform5x2,performance_ci)))+
  scale_color_manual(name="Kernel",values=c("firebrick1","sienna2","goldenrod2","plum3","skyblue3"))

 
plot_grid(p,q,ncol=2,rel_widths = c(1,1.2))

##### Comprovar que la predicció és bona en tota el domini
for(j in 1:5) {
  par(mfrow=c(4,4))

  for(drug in DRUGS){
    print(plot( DATA[train_idx,drug],FIT[train_idx,drug,j], cex=0.5,pch=20,main=paste(drug, "training: actual vs fitted"),
          # xlab = "Actual",ylab="Fitted") ## Fitted values
          xlab = "Actual",ylab="Fitted"),fill="steelblue3") ## Fitted values
    print(text(x = 0, y = 2, # Coordinates
               label =  expression("R" ^ 2 * "=") ))
    print(text(x = 0.5, y = 1.95, # Coordinates
               label =  round( cor( DATA[train_idx,drug] ,FIT[train_idx,drug,j],
                                    use = "pairwise.complete.obs")^2,digits=2 )))
    print(plot( DATA[te_idx,drug],PREDS[te_idx,drug,j],  cex=0.5,pch=20, main=paste(drug, "test: actual vs predicted"),
          xlab = "Actual",ylab="Predicted"),fill="steelblue3") ## Predicted values
    if(drug=="TPV") {
      print(text(x = 2, y = 0, # Coordinates
                 label =  expression("R" ^ 2 * "=") ))
      print(text(x = 2.5, y = -0.01, # Coordinates
                 label =  round( cor(   DATA[te_idx,drug],PREDS[te_idx,drug,j],
                                        use = "pairwise.complete.obs")^2,digits=2 )))

    } else {
      print(text(x = 0, y = 2, # Coordinates
                 label =  expression("R" ^ 2 * "=") ))
      print(text(x = 0.5, y = 1.95, # Coordinates
                 label =  round( cor(   DATA[te_idx,drug],PREDS[te_idx,drug,j],
                                        use = "pairwise.complete.obs")^2,digits=2 )))
    }

  }
}
par(mfrow=c(1,1))

### Performance (test alternatiu)
alt_te_idx <-  which(Partition ==  "Additional test")

performance_alt <- matrix(0,nrow=length(DRUGS),ncol=6)
rownames(performance_alt) <- DRUGS
colnames(performance_alt)  <- c("Spectrum1","Spectrum2","Spectrum3","Intersect","Interactions","Nte")

for(drug in DRUGS) {
  print(drug)
  Ytest <- DATA[alt_te_idx,drug]
  cc <- complete.cases(Ytest) &  !is.infinite(Ytest)
  performance_alt[drug,"Nte"] <- sum(cc)
  for(j in 1:5)    performance_alt[drug,j] <- nmse(Ytest[cc],PREDS[alt_te_idx[cc],drug,j])
}

### amb CIs. Comprovar que la predicció és bona en tot el domini
performance_ci_alt <- array(0,dim=c(length(DRUGS),3,5),
                        dimnames =  list(DRUGS,
                                         c("Mean","2.5%","97.5%"),
                                         c("Spectrum1","Spectrum2","Spectrum3","Intersect","Interactions")))
for(drug in DRUGS) {
  print(drug)
  Ytest <- DATA[alt_te_idx,drug]
  cc <- complete.cases(Ytest) &  !is.infinite(Ytest)
  for(j in 1:5)    performance_ci_alt[drug,,j] <- Boots_CI(target=Ytest[cc], pred=PREDS[alt_te_idx[cc],drug,j],
                                                       index = "nmse", nboots=30000, confidence = 95)
}

#### Taula
PERF_alt <- performance_alt
for(j in 1:5) PERF_alt[,j] <- apply(round(performance_ci_alt[,,j],digits=3),1,function(x)paste0(x[1]," (",x[2],", ",x[3],")"))

write.csv(PERF_alt,file = "Results/performance_altest.csv")




### Importàncies
imps_first20 <- array(0,dim=c(length(DRUGS),20,5),
      dimnames =  list(DRUGS,
                       1:20,
                       c("Spectrum1","Spectrum2","Spectrum3","Intersect","Interactions")))

cumsums <- performance

for(drug in DRUGS) {
  A <- plotImp(IMPS$Intersect[  ,drug],nfeat=20,relative = FALSE,
                                   absolute=FALSE,leftmargin=6,main=paste("Intersect kernel",drug))
  B <- plotImp(IMPS$Spectrum1[ ,drug], relative = FALSE,absolute=FALSE,leftmargin=6,main="Spectrum1 kernel")
  C <- plotImp(IMPS$Spectrum2[ ,drug],nfeat=20,relative = TRUE,
                                   absolute=FALSE,leftmargin=6,main=paste("Spectrum2 kernel",drug))
  D <- plotImp(IMPS$Spectrum3[ ,drug],nfeat=20,relative = TRUE,
                                   absolute=FALSE,leftmargin=6,main=paste("Spectrum3 kernel",drug))
  E  <- plotImp(IMPS$Interactions[ ,drug],nfeat=20,relative = TRUE,absolute=TRUE,leftmargin=6,main="Interactions1 kernel")

  imps_first20[drug,,4] <- A$first_features
  imps_first20[drug,,1] <- B$first_features
  imps_first20[drug,,2] <- C$first_features
  imps_first20[drug,,3] <- D$first_features
  imps_first20[drug,,5] <- E$first_features

  cumsums[drug,4] <- A$cumsum/sum(abs(IMPS$Intersect[ ,drug]))
  cumsums[drug,1] <- B$cumsum/sum(abs(IMPS$Spectrum1[ ,drug]))
  cumsums[drug,2] <- C$cumsum/sum(abs(IMPS$Spectrum2[ ,drug]))
  cumsums[drug,3] <- D$cumsum/sum(abs(IMPS$Spectrum3[ ,drug]))
  cumsums[drug,5] <- E$cumsum/sum(abs(IMPS$Interactions[ ,drug]))
}

### intersect
round(cumsums,digits=3)
mutations_intersect_20 <- matrix("",nrow=length(PrPositions),ncol=length(DRUGS))
rownames(mutations_intersect_20) <- PrPositions
colnames(mutations_intersect_20) <- DRUGS

for(drug in DRUGS) {
  pos <- sub('_.*', '', imps_first20[drug,,4])
  alleles <-  sub('.*_', '', imps_first20[drug,,4])
  for(al in 1:20)   mutations_intersect_20[pos[al],drug] <- paste0(mutations_intersect_20[pos[al],drug],alleles[al])

}

## On són les posicions 46 i 71?
for(drug in DRUGS)print(grep("P46_",names(IMPS$Intersect[order(abs(IMPS$Intersect[,drug]),decreasing=T),drug])))
for(drug in DRUGS)print(grep("P71_",names(IMPS$Intersect[order(abs(IMPS$Intersect[,drug]),decreasing=T),drug])))


mutations_intersect_20 <- apply(as.data.frame(mutations_intersect_20),c(1,2),strSort)
rownames(mutations_intersect_20) <- paste0(sub("^P", "", rownames(mutations_intersect_20)),WTSeq)
mutations_intersect_20 <- mutations_intersect_20[rowSums(mutations_intersect_20=="")<length(DRUGS),]
mutations_intersect_20


## Interactions
IMPS_summed <- list()
IMPS_summed$Interactions <- aggregate_imp(abs(IMPS$Interactions^2),
                                          lev= unique(substr(rownames(IMPS$Interactions),1,nchar(rownames(IMPS$Interactions))-3)),samples="cols")

mutations_interactions_20 <- matrix(0,nrow=length(unique(as.vector(imps_first20[,,5] ))),ncol=length(DRUGS))
rownames(mutations_interactions_20) <- unique(as.vector(imps_first20[,,5] ))
colnames(mutations_interactions_20) <- DRUGS
for(drug in DRUGS)   mutations_interactions_20[imps_first20[drug,,5],drug] <- 1
sort(rowSums(mutations_interactions_20),decreasing = T)[1:20]
mutations_interactions_20 <- as.data.frame(mutations_interactions_20)
mutations_interactions_20$pos1 <- "Neutral"
mutations_interactions_20$pos2 <- "Neutral"

mutations <-   list(major= paste0("P",c(23, 24, 30, 32, 33, 46, 47, 48, 50, 53, 54, 73, 76, 82, 84, 88, 90)),
                    minor = paste0("P",c(10, 11, 13, 20, 34, 35, 36, 43, 45, 55, 58, 60, 63, 71, 74, 75, 77, 79, 83, 85, 89, 91, 93, 95)))
for(posxpos in rownames(mutations_interactions_20)) {
  inmajor <- unlist(strsplit( (posxpos),split="_")) %in% mutations$major
  inminor <- unlist(strsplit( (posxpos),split="_")) %in% mutations$minor
  mutations_interactions_20[posxpos,9:10][,inmajor] <-"Major"
  mutations_interactions_20[posxpos,9:10][,inminor] <-"Minor"
}
sort(rowSums(mutations_interactions_20[,DRUGS]),decreasing = T)


DATAdum_delete <- rearrange_kint(Kint$feat_space,rows=rownames(DATA),AA=AA,PrPositions = PrPositions)

colors_wt <- colourx(IMPS$Interactions,WTSeq,DATAdum_delete)

mutations_interactions_20_alleles <- mutations_interactions_20
mutations_interactions_20_alleles[,DRUGS] <- ""
for(drug in DRUGS) {
  posxpos <- rownames(mutations_interactions_20)[mutations_interactions_20[,drug]>0]
  for(pos in posxpos) {
        searchfor <- grep(pos,rownames(IMPS$Interactions),value = T)
        if(length(searchfor)<2) stop("error")
        # allelesmore <-  names(sort( (IMPS$Interactions[searchfor,drug]),decreasing = T)[1:2])
        allelesmore <-  names(sort( (IMPS$Interactions[searchfor,drug]),decreasing = T)[1:2])
        alleless <-  names(sort( (IMPS$Interactions[searchfor,drug]),decreasing = F)[1:2])
        # alleles <-  names(sort(abs(IMPS$Interactions[searchfor,drug]),decreasing = T)[1:5])
        allelesmore <-         unique(substr(allelesmore, nchar(allelesmore)-1,nchar(allelesmore)))
        alleless <-         unique(substr(alleless, nchar(alleless)-1,nchar(alleless)))
        mutations_interactions_20_alleles[pos,drug] <- paste(paste(alleless,collapse = "/"),paste(allelesmore,collapse = "/"),sep=";")
  }
}

changernames <-  rownames(mutations_interactions_20_alleles)
changernames <- unlist(strsplit(changernames,"_",))
# cosa <- which(II[,first20,drug]>threshold,arr.ind=T)
# cosacosa <- II[,,drug]
# cosacosa[] <- ""
changernames <- data.frame(row=sub("P","",changernames[seq(1,length(changernames),2)]),col=  sub("P","",changernames[seq(2,length(changernames),2)]))
changernames <- paste(paste0(WTSeq[as.numeric(changernames[,1])],changernames[,1]),paste0(WTSeq[as.numeric(changernames[,2])],changernames[,2]),sep="_")
rownames(mutations_interactions_20_alleles) <- changernames
mutations_interactions_20_alleles[order(rownames(mutations_interactions_20)),]



## Spectrum-3
COINC_DRUGS_ALT <-   (IMPS$Spectrum3) [sort(unique(as.vector(imps_first20[,,3]))) ,]

for(drug in DRUGS) {
  whichare <- rownames(COINC_DRUGS_ALT) %in% imps_first20[drug ,,3]
  COINC_DRUGS_ALT[ !whichare,drug] <- 0
}

heatmapcols <- rep("wheat2",nrow(COINC_DRUGS_ALT))
heatmapcols[which(rownames(COINC_DRUGS_ALT) %in% paste(names(WTSeq),WTSeq,sep="_"))] <- "skyblue"
heatmapcols[which( vect_grep(rownames(COINC_DRUGS_ALT),WTSeqseq))] <- "skyblue"

heatmap(  (COINC_DRUGS_ALT),scale = "column",main="Top 20 3-mers",  
          Colv =NA,
        RowSideColors = heatmapcols)
  dev.off()


  
  ### Comparació de rankings (feature importances)

Kpi4 <- Kendall( IMPS$Intersect , samples.in.rows = FALSE)

kint_kendall <- kPCA(Kpi4 ,center=T, plot=1:2, name_leg = "SVM error",
                     y=performance[,"Intersect"],pos_leg = "none",
                    colors = c("#66CDAA","orange","#CD3278"),title="Intersect SVM Importances")

 
#### Sp3

Kpi3 <- Kendall(IMPS$Spectrum3 , samples.in.rows = FALSE)
Kpi3
 
sp3_kendall <- kPCA(Kpi3,center=T ,plot=1:2,
                    name_leg = "SVM error",
                    y=performance[,"Spectrum3"],  colors = c("#66CDAA","orange","#CD3278"),title="Spectrum k=3 SVM Importances")

 

#### Interaction
Kpi5 <- Kendall(IMPS$Interactions,samples.in.rows = FALSE)
 
mixlong1_kendall <- kPCA(Kpi5,center=T,
                         plot=1:2, name_leg = "SVM error",
                         y=performance[,"Interactions"],pos_leg = "none",
                         colors = c("#66CDAA","orange","#CD3278"),title="Interaction SVM Importances")
mixlong1_kendall$plot

legend <- ggplotGrob(sp3_kendall$plot + theme(legend.position = "bottom"))$grobs
legend <- legend[[which(sapply(legend, function(x) x$name) == "guide-box")]]

plot_grid(
  plot_grid(
    plot_grid(kint_kendall$plot +geom_text(label=DRUGS,hjust= "inward", vjust= "inward"),sp3_kendall$plot+geom_text(label=DRUGS,hjust= "inward", vjust= "inward")+ ggplot2::theme(legend.position="none"), nrow = 1, ncol = 2, labels=c("A","B")),
    plot_grid(NULL,   mixlong1_kendall$plot+geom_text(label=DRUGS,hjust= "inward", vjust= "inward"), NULL, nrow = 1, rel_widths = c(0.5, 1, 0.5),labels = "C"),
    nrow = 2
  ), legend,nrow=2,rel_heights = c(0.9,0.1))


