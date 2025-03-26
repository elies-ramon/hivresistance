############################
##### KERNEL MATRICES ######
############################

# Llibreries #
library(kerntools)


# Càrrega de les funcions #
source("helper_functions.R")


# Càrrega dels indexes #
load("Partitions/PI_IDX")


# Càrrega de les posicions de resistència #
WTSeq <- read.table("Data/outputPI",
                    header=FALSE, col.names=paste0("P",1:99),sep="",colClasses = c("character"),
                    stringsAsFactors=TRUE)
WTSeqseq <- toSeq(WTSeq)


# Càrrega de dades #

## Sense mescles d'aa
load("Data/PI_Seq_clean") # dades
load("LD/PI_DIST.dat") # distancies
diag(DIST99) <- NA


for (i in 1:30)   {

  DATA.l <-   DATA.list[[i]]

  # Particionat
  DATA <- droplevels(DATA.l[c(IDX$work,IDX$additional_test_subtype),])
  dim(DATA)  #2212  115
  DATAseq <- toSeq(DATA[,PrPositions])
  DATAseq <- cbind(DATA[,setdiff(colnames(DATA),PrPositions)],DATAseq)

  # Càlcul del kernel #

  ## Spectrum kernel
  SP1 <- Spectrum(DATAseq$DATAseq, alphabet=AA, l=1, feat_space=TRUE)
  SP2 <- Spectrum(DATAseq$DATAseq, alphabet=AA, l=2, feat_space=TRUE)
  SP3 <- Spectrum(DATAseq$DATAseq, alphabet=AA, l=3, feat_space=TRUE)

  save(SP1,file=paste0("Kernel_matrices/PI_SP1_",i,".dat"))
  save(SP2,file=paste0("Kernel_matrices/PI_SP2_",i,".dat"))
  save(SP3,file=paste0("Kernel_matrices/PI_SP3_",i,".dat"))

  dim(SP1$feat_space)
  dim(SP2$feat_space)
  dim(SP3$feat_space)

  SP2$feat_space <- desparsify(SP2$feat_space,2)
  SP3$feat_space <- desparsify(SP3$feat_space,2)
  dim(SP2$feat_space)
  dim(SP3$feat_space)
  
  ## Dirac kernel
  Kov <- Dirac (DATA[,PrPositions],comp = "sum" ,feat_space=TRUE )
  save(Kov,file=paste0("Kernel_matrices/PI_Kov_",i,".dat"))

  dim(Kov$K)
  dim(Kov$feat_space)
  
  load(file=paste0("Kernel_matrices/PI_Kov_",i,".dat"))

  ## Interactions-Dirac kernel (parelles d'aa)
  #### Refer: Polato et al 2018 ("monotone conjuctive kernel")
  Long1 <- choose(Kov$K,2)
  DATAlong1 <- apply(apply(DATA[,PrPositions],2,as.character),1,pos_interac)
  DATAlong1 <- t(DATAlong1)
  cols <- t(combn(PrPositions,2))
  colnames(DATAlong1) <- paste(cols[,1],cols[,2],sep="_")
  DATAlong1 <- dummy_data(DATAlong1)
  LONG1 <- list(K=Long1,feat_space=DATAlong1)

  save(LONG1,file=paste0("Kernel_matrices/PI_LONG1_",i,".dat"))

  # Retallant interacions <=15a
  load(file=paste0("Kernel_matrices/PI_LONG1_",i,".dat"))
  pos15 <- which(DIST99>15 | is.na(DIST99),arr.ind = T)
  todelete <- c()
  todelete <- c()
  for(j in 1:nrow(pos15)) {
    if(pos15[j,1]<pos15[j,2]) {
      todelete <- append(todelete,paste0(PrPositions[pos15[j,1]],"_",PrPositions[pos15[j,2]],"_"))
    } else {
      todelete <- append(todelete,paste0(PrPositions[pos15[j,2]],"_",PrPositions[pos15[j,1]],"_"))
    }
  }

  for(posx in todelete) {
    LONG1$feat_space[,grep(posx,colnames(LONG1$feat_space))] <- 0
  }
  LONG1$feat_space <- desparsify(LONG1$feat_space)
  dim(LONG1$feat_space )
  LONG1$K  <- Linear(LONG1$feat_space)

  save(LONG1,file=paste0("Kernel_matrices/PI_LONG_15a",i,".dat"))

  gc()
}


## Mitjanes

N <- length(c(IDX$work,IDX$additional_test_subtype))

Ksp1m <- Ksp2m <- Ksp3m <- Kovm <- Kov2m <- Long1m <- Long2m <- array(0, dim=c(N,N,30))


for(i in 1:30) {
  load(file = paste0("Kernel_matrices/PI_SP1_",i,".dat"))
  Ksp1m[,,i] <- SP1$K
  rm(SP1)
  gc()
}

Ksp1m <- rowMeans(Ksp1m,dims = 2)
save(Ksp1m,file= "Kernel_matrices/PI_SP1_mean.dat")
rm(Ksp1m)


for(i in 1:30) {
  load(file = paste0("Kernel_matrices/PI_SP2_",i,".dat"))
  Ksp2m[,,i] <- SP2$K
  rm(SP2)
  gc()
}

Ksp2m <- rowMeans(Ksp2m,dims = 2)
save(Ksp2m,file= "Kernel_matrices/PI_SP2_mean.dat")
rm(Ksp2m)


for(i in 1:30) {
  load(file = paste0("Kernel_matrices/PI_SP3_",i,".dat"))
  Ksp3m[,,i] <- SP3$K
  rm(SP3)
  gc()
}

Ksp3m <- rowMeans(Ksp3m,dims = 2)
save(Ksp3m,file= "Kernel_matrices/PI_SP3_mean.dat")
rm(Ksp3m)


for(i in 1:30) {
  load(file = paste0("Kernel_matrices/PI_Kov_",i,".dat"))
  Kovm[,,i] <- Kov$K
  rm( Kov)
  gc()
}

Kovm <- rowMeans(Kovm,dims = 2)
save(Kovm,file= "Kernel_matrices/PI_Kov_mean.dat")
rm(Kovm)


for(i in 1:30) {
  load(file = paste0("Kernel_matrices/PI_LONG_15a",i,".dat"))
  Long1m[,,i] <- LONG1$K
  rm(LONG1)
  gc()
}

Long1m <- rowMeans(Long1m,dims = 2)
save(Long1m,file= "Kernel_matrices/PI_LONG_15a_mean.dat")
rm(Long1m)


###################

## Amb mescles d'aa

# Càrrega dels indexes #
load("Partitions/PI_IDX")


# Càrrega de dades #
load("Data/PI_MultiA_clean")

# Particionat
DATA.train <- DATA[c(IDX$work,IDX$additional_test_subtype),]
dim(DATA.train)   

te_idx <- DATA[IDX$test,"SeqID"]
te_add_subtype <-  DATA[IDX$additional_test_subtype,"SeqID"]
anyDuplicated(te_idx)
anyDuplicated(te_add_subtype)


# Càlcul del kernel #

## Intersect kernel
Kint <- Intersect (DATA.train[,PrPositions],elements = AA,comp = "mean" ,feat_space=TRUE ) # 10GB de RAM. No és viable.
dim(Kint$K)
dim(Kint$feat_space)
#  2212 2212
# [1] 2212   20   99

save(Kint,  file= "Kernel_matrices/PI_Kint.dat")


## Spectrum-mixtures
load("Kernel_matrices/PI_SP1_1.dat")
load("Kernel_matrices/PI_SP2_1.dat")
load("Kernel_matrices/PI_SP3_1.dat")

DATAseqfs <- SP3$feat_space
DATAgenerate <- apply(DATA.train[,PrPositions],c(1,2),nchar)
 sum(apply(DATAgenerate,1,prod))
 
DATAforseq <- apply(DATA.train[,PrPositions],1:2,as.character)
L <- 1
idx <-  which(apply(DATAgenerate,1,prod)>1)
alphabet <- permute_rep(alphabet = AA, l = L)

  for(i in idx) {
    minidata <- matrix(substring(DATAforseq[i, ],1,1),nrow= apply(DATAgenerate,1,prod)[i],ncol=99,byrow = T)
    chars <-  matrix(c(substring(DATAforseq[i, ],1,1),substring(DATAforseq[i, ],2,2),substring(DATAforseq[i, ],3,3),
                     substring(DATAforseq[i, ],4,4)),
                     nrow= 4,ncol=99,byrow = T)
    colnames(minidata) <-  colnames(chars) <- (PrPositions)
    combs <- DATAgenerate[i,which(DATAgenerate[i,]>1),drop=F]
    cols <- colnames(combs)
    combs <- expand.grid( lapply(combs,function(x)1:x))
    colnames(combs) <- cols
    for(j in 1:ncol(combs)) combs[,j] <-  chars[combs[,j],colnames(combs)[j]]
    minidata[,colnames(combs)] <- as.matrix(combs)
    rm(chars,combs)
    minidata <- toSeq(minidata)
    count_table <- searchSubs(minidata, fixed = alphabet, overlap = TRUE)
    DATAseqfs[i,] <- colMeans(count_table)
    rm(minidata,count_table)
    gc()
  }

dim(DATAseqfs)
mixSP3 <- SP3
mixSP3$K <- Linear(DATAseqfs)
mixSP3$feat_space <- DATAseqfs
save(mixSP3,file="Kernel_matrices/PI_mixSP3.dat")


## Interactions-Intersect kernel (parelles d'aa)
Mixlong1 <- choose(length(PrPositions)*Kint$K,2)
DATAdum_delete <- matrix(Kint$feat_space, nrow(Kint$K), 99*20)
rownames(DATAdum_delete) <- rownames(Kint$K)
colnames(DATAdum_delete) <- paste(rep(paste0("P",1:99),each=20),rep(AA,99),sep="_")
DATAdum_delete <- desparsify(DATAdum_delete,2)
dim(DATAdum_delete)
DATAlong1 <- matrix(0,nrow=nrow(DATAdum_delete),ncol=choose(ncol(DATAdum_delete),2))
for(i in 1:nrow(DATAdum_delete))   DATAlong1[i,] <- combn(DATAdum_delete[i,],2,FUN=prod)
cols <- t(combn(colnames(DATAdum_delete),2))
cols1 <-  sub(".*_", "", cols)
cols2 <-  sub("_.*", "", cols)
colnames(DATAlong1) <- paste(paste(cols2[,1],cols2[,2],sep="_"), paste0(cols1[,1],cols1[,2]) ,sep="_")
Mixlong1 <- list(K=Mixlong1,feat_space=DATAlong1)
save(Mixlong1,file="Kernel_matrices/PI_MIXLONG1.dat")

pos15 <- which(DIST99>15 | is.na(DIST99),arr.ind = T)
todelete <- c()
for(i in 1:nrow(pos15)) {
  if(pos15[i,1]<pos15[i,2]) {
    todelete <- append(todelete,paste0(PrPositions[pos15[i,1]],"_",PrPositions[pos15[i,2]],"_"))
  } else {
    todelete <- append(todelete,paste0(PrPositions[pos15[i,2]],"_",PrPositions[pos15[i,1]],"_"))
  }
}

for(posx in todelete) {
  Mixlong1$feat_space[,grep(posx,colnames(Mixlong1$feat_space))] <- 0
}
Mixlong1$feat_space <- desparsify(Mixlong1$feat_space)
dim(Mixlong1$feat_space )
Mixlong1$K  <- Linear(Mixlong1$feat_space)
save(Mixlong1,file="Kernel_matrices/PI_MIXLONG_15a.dat")