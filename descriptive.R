########################
##### DESCRIPTIVA ######
########################

# Declaració de funcions #

library(kerntools)
library(maotai)
library(ggplot2)
library(viridis)
library(cowplot)
library(reshape2)


source("helper_functions.R")

# Càrrega dels indexes #
load("Partitions/PI_IDX")

# Càrrega de dades #
load("Data/PI_MultiA_clean")
rownames(DATA) <- 1:nrow(DATA)



# N #

## Tr,te & test alternatiu

te_add_subtype <-  DATA[IDX$additional_test_subtype,"SeqID"]
te_idx <- DATA[ IDX$test,"SeqID"]

## Tr per fàrmac
for(drug in DRUGS) cat(paste0(drug,": " ,sum(complete.cases(DATA[as.numeric(IDX$work[!(IDX$work %in% IDX$test)]),drug])),"\n"))


## Total
DATA <- droplevels(DATA[as.numeric(c(IDX$work,IDX$additional_test_subtype)),])
dim(DATA)  #2212  115
nrow(DATA)-length(IDX$additional_test_subtype) -length(IDX$test) # Training 1680
length(IDX$additional_test_subtype) # Test alternatiu: 125
length(IDX$test) # Test: 407

## Resistència per fàrmac ##
colors <- plasma(8)
plotresistance <- data.frame(DATA[,DRUGS])
min(plotresistance,na.rm = T)
max(plotresistance,na.rm = T)

A <- ggplot(plotresistance, aes(x=FPV)) +  geom_density(color=colors[1],fill=colors[1],alpha=0.4) +  theme_bw() + xlim(-1,3) + ggtitle("Log-Resistance distribution")
B <- ggplot(plotresistance, aes(x=ATV)) +  geom_density(color=colors[2],fill=colors[2],alpha=0.4) +  theme_bw() + xlim(-1,3)
C <- ggplot(plotresistance, aes(x=IDV)) +  geom_density(color=colors[3],fill=colors[3],alpha=0.4) +  theme_bw() + xlim(-1,3)
D <- ggplot(plotresistance, aes(x=LPV)) +  geom_density(color=colors[4],fill=colors[4],alpha=0.4) +  theme_bw() + xlim(-1,3)
E <- ggplot(plotresistance, aes(x=NFV)) +  geom_density(color=colors[5],fill=colors[5],alpha=0.4) +  theme_bw() + xlim(-1,3)
FF <- ggplot(plotresistance, aes(x=SQV)) +  geom_density(color=colors[6],fill=colors[6],alpha=0.4) +  theme_bw() + xlim(-1,3)
G <- ggplot(plotresistance, aes(x=TPV)) +  geom_density(color=colors[7],fill=colors[7],alpha=0.4) +  theme_bw() + xlim(-1,3)
H <- ggplot(plotresistance, aes(x=DRV)) +  geom_density(color=colors[8],fill=colors[8],alpha=0.4) +  theme_bw() + xlim(-1,3)

plot_grid(A,B,C,D,E,FF,G,H,nrow=8,rel_heights = c(1.2,rep(1,7)))

## Mutacions per sequencia ##

mutsperseq <- Intersect(rbind(WTSeq , DATA[,PrPositions]),comp = "mean",elements = AA)
mutsperseq <- 99*mutsperseq
rownames(mutsperseq) <- colnames(mutsperseq) <- c("wt",rownames(DATA))
mutsperseq <- (99-mutsperseq[1,-1])
hist(mutsperseq,main="Number of mutations respect to the WT sequence",xlab="N mutations")
summary(mutsperseq[DATA$Subtype=="B"])
summary(mutsperseq[DATA$Subtype!="B"])


# Mutacions x fàrmac #
mutsperdrug <- data.frame( mutsperseq,DATA[, DRUGS])
colnames(mutsperdrug) <- c("Mutations",DRUGS)
mutsperdrug <- data.frame(Mutations=rep(mutsperdrug$Mutations),melt(mutsperdrug[,DRUGS]))
colnames(mutsperdrug)[2] <- c("Drug")

ggplot(mutsperdrug, aes(x = Mutations, color = Drug)) +geom_point(aes(y = value),size=0.8) + ggtitle("Log-resistance vs number of mutations") +
  scale_color_manual(values = colors)+facet_wrap(~Drug,  ncol=2,) + theme_bw() + ylab("Log-resistance") +xlab("N mutations") +
  theme(legend.position = "none",strip.background = element_rect(colour="white", fill="white"))


# Mescles d'aminoàcids #
multi_by_sample <- rowSums(apply(DATA[,PrPositions],c(1,2),nchar)>1)
sum(multi_by_sample>0)/nrow(DATA) #  0.5831826 al menys tenen 1 posició amb barreja d'aa

par(mfrow=c(1,2))

hist(multi_by_sample,main="Allelic mixtures per sample",xlab = "Number of affected positions",
     ylab = "Absolute frequency", border="white",col="#7e5dd9")
summary(multi_by_sample)
multi_by_position <- colSums(apply(DATA[,PrPositions],c(1,2),nchar)>1)
sort(multi_by_position,decreasing = T)[1:20]
summary( (multi_by_position  ))

# Frequències al·lèliques#
sum(sapply(DATA[,PrPositions],nlevels)>1)/99  # 0.979798 són polimòrfiques. Excepcions: P1_P i P42_W

load("LD/Freqs_mixt.dat") # amb mescles d'aa

boxplot(sort(colMeans(FREQ0_mixt),decreasing = T)[1:99],main="Major allele frequency",
        ylab="Relative Frequency",pch=16,border="#7e5dd9",col="white")

sort(sort(colMeans(FREQ0_mixt),decreasing = T)[1:99])
  
# Kernel Two-sample Test with Maximum Mean Discrepancy
Partition <- rep("Training",nrow(DATA))
Partition[which( DATA$SeqID %in% te_idx)] <- "Test"
Partition[which( DATA$SeqID %in% te_add_subtype)] <- "Additional test"

## Additional test vs test
load("Kernel_matrices/PI_Kint.dat")
idx <- which(Partition!="Training")
mmd2test(Kint$K[idx,idx],label = Partition[idx])

Partition[which( DATA$SeqID %in% te_add_subtype)] <- "Additional test"

## Training vs test
idx <- which(Partition!="Additional test")
mmd2test(Kint$K[idx,idx],label = Partition[idx])
 
# Training vs Additional test
idx <- which(Partition!="Test")
mmd2test(Kint$K[idx,idx],label = Partition[idx])
