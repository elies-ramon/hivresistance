########################
##### kernel PCAs ######
########################


# Llibreries #
library(kerntools)


# Declaració de funcions i constants #
source("helper_functions.R")



# Càrrega dels indexes #
load("Partitions/PI_IDX")



# Càrrega de dades #
load("Data/PI_Seq_clean") # dades
load("LD/PI_LD.dat") # desequilibri lligament


PCS_S1 <- PCS_S2 <- PCS_S3 <- PCS_OV <- list()

## Sense mescles d'aa
for (i in 1:30)   {

  DATA.l <-   DATA.list[[i]]

  # Particionat #
  DATA <- droplevels(DATA.l[c(IDX$work,IDX$additional_test_subtype),])
  dim(DATA)  #2212  115
  DATAseq <- toSeq(DATA[,PrPositions])
  DATAseq <- cbind(DATA[,setdiff(colnames(DATA),PrPositions)],DATAseq)

  te_idx <- DATA.l[IDX$test,"SeqID"]
  te_add_subtype <-  DATA.l[IDX$additional_test_subtype,"SeqID"]
  anyDuplicated(te_idx)
  anyDuplicated(te_add_subtype)

  Partition <- rep("Training",nrow(DATAseq))
  Partition[which( DATAseq$SeqID %in% te_idx)] <- "Test"
  Partition[which( DATAseq$SeqID %in% te_add_subtype)] <- "Additional test"

  # Càrrega matrius #
  load(paste0("Kernel_matrices/PI_SP1_",i,".dat"))
  load(paste0("Kernel_matrices/PI_SP2_",i,".dat"))
  load(paste0("Kernel_matrices/PI_SP3_",i,".dat"))
  load(paste0("Kernel_matrices/PI_Kov_",i,".dat"))
  load(paste0("Kernel_matrices/PI_LONG1_",i,".dat"))

  # Estudi diversitat #
  #### En aquesta part es pot afegir els grups != de B per vore en cauen de sa PCA.

  ## Spectrum kernel
  ### abund kmers
  sp1_abund <- colMeans(SP1$feat_space)
  sp2_abund <- colMeans(SP2$feat_space)
  sp3_abund <- colMeans(SP3$feat_space)

  SP2$feat_space <- desparsify(SP2$feat_space,2)
  SP3$feat_space <- desparsify(SP3$feat_space,2)
  dim(SP2$feat_space)
  # 2212  346
  dim(SP3$feat_space)
  # 2212 1797

  SP1$feat_space_ce <- scale(SP1$feat_space,center=T,scale=F)
  SP2$feat_space_ce <- scale(SP2$feat_space,center=T,scale=F)
  SP3$feat_space_ce <- scale(SP3$feat_space,center=T,scale=F)

  ### Càlcul kernel centrat
  Ksp1 <- Linear(SP1$feat_space_ce)
  Ksp2 <- Linear(SP2$feat_space_ce)
  Ksp3 <- Linear(SP3$feat_space_ce)

  ### String kernel PCA
  pca_s1 <- kPCA(Ksp1,center=F,plot=1:2,colors = "grey99")
  pca_s2 <- kPCA(Ksp2,center=F,plot=1:2,colors = "grey99" )
  pca_s3 <- kPCA(Ksp3,center=F,plot=1:2,colors = "grey99" )
  # pca_s4 <- kPCA(Ksp4,center=F,plot=1:2,colors = "grey99" )

  #### Contribucions de cada kmer als PCs.
  pcs_s1 <- kPCA_imp(SP1$feat_space_ce,center = FALSE,projected=pca_s1$projection,secure=FALSE)
  pcs_s2 <- kPCA_imp(SP2$feat_space_ce,center = FALSE,projected=pca_s2$projection,secure=FALSE)
  pcs_s3 <- kPCA_imp(SP3$feat_space_ce,center = FALSE,projected=pca_s3$projection,secure=FALSE)

  PCS_S1[[i]] <- pcs_s1$loadings[1:2,]
  PCS_S2[[i]] <- pcs_s2$loadings[1:2,]
  PCS_S3[[i]] <- pcs_s3$loadings[1:3,]


  pc1_list2 <- plotImp(PCS_S2[[i]][1,], absolute=FALSE, relative = FALSE,  main="PC1",
                        leftmargin = 5,  color="pink",nfeat=20)
  pc2_list2 <- plotImp(PCS_S2[[i]][2,], absolute=FALSE, relative = FALSE,  main="PC1",
                        leftmargin = 5,  color="pink",nfeat=20)
  pc1_list3 <- plotImp(PCS_S3[[i]][1,], absolute=FALSE, relative = FALSE,  main="PC1",
                        leftmargin = 5,  color="pink",nfeat=20)

  pc2_list3 <- plotImp(PCS_S3[[i]][2,], absolute=FALSE, relative = FALSE,  main="PC1",
                       leftmargin = 5,  color="pink",nfeat=20)

 

  ## Dirac kernel
  ### abund
  dummies_abund <- colMeans(Kov$feat_space) # en tant per 1

  ### Plantejar-se si esborrar les posicions constants
  # DATAdum <- DATAdum[,-which(dummies_abund==1)] # Ho he fet en el passat però em sembla una liada, i més comparant amb string kernels

  ### Càlcul kernel centrat
  DATAdum_ce <- scale(Kov$feat_space,center=TRUE,scale=FALSE)
  Kov <- Linear(DATAdum_ce)

  ### Dirac kernel PCA
  pca_ov <- kPCA(Kov,center=FALSE,plot=1:2,colors = "grey99")


  #### Contribucions de cada al·lel als PCs
  pcs_ov <- kPCA_imp(DATAdum_ce,center = FALSE,projected=pca_ov$projection,secure=FALSE)
  PCS_OV[[i]] <- pcs_ov$loadings[1:2,]

  pc1_list <- plotImp(pcs_ov$loadings[1,], absolute=FALSE, relative = FALSE,  main="PC1",
                        leftmargin = 5,  color="pink",nfeat=20)
  pc2_list <- plotImp(pcs_ov$loadings[2,], absolute=FALSE, relative = FALSE,  main="PC1",
                        leftmargin = 5,  color="pink",nfeat=20)

 
  
    # ## Interactions-Dirac kernel
  pca_long1 <- kPCA(Long1,center=TRUE,plot=1:2,colors = "grey99")

}

# mitjana de 30 mostres
## Dirac kernel
load("Kernel_matrices/PI_Kov_mean.dat")

### Dirac kernel PCA
pca_ov <- kPCA(Kovm,center=TRUE,plot=1:2,colors = "grey99")


#### Contribucions de cada al·lel als PCs
pcs_ov <- kPCA_imp(DATAdum_ce,center = FALSE,projected=pca_ov$projection,secure=FALSE)
PCS_OV[[i]] <- pcs_ov$loadings[1:2,]

pc1_list <- plotImp(pcs_ov$loadings[1,], absolute=FALSE, relative = FALSE,  main="PC1",
                    leftmargin = 5,  color="pink",nfeat=20)
pc2_list <- plotImp(pcs_ov$loadings[2,], absolute=FALSE, relative = FALSE,  main="PC1",
                    leftmargin = 5,  color="pink",nfeat=20)


 
## Amb mescles d'aa

# Càrrega de dades #

load("Data/PI_MultiA_clean")
load("Kernel_matrices/PI_Kint.dat")

# Particionat
rownames(DATA) <- 1:nrow(DATA)
DATA.train <- droplevels(DATA[c(IDX$work,IDX$additional_test_subtype),])
dim(DATA.train)  #2212  115
summary(factor(DATA.train$Type))

te_idx <- DATA[IDX$test,"SeqID"]
te_add_subtype <-  DATA[IDX$additional_test_subtype,"SeqID"]
anyDuplicated(te_idx)
anyDuplicated(te_add_subtype)


## Intersect kernel
### abund kmers
dummies_abund <- colMeans(Kint$feat_space)


#Reorganitzar + Esborrar posicions que apareixen 0 voltes:
DATAdum_delete <- rearrange_kint(Kint$feat_space)


### Càlcul kernel centrat
DATAdum_ce <- scale(DATAdum_delete,center=TRUE,scale=FALSE)
Kint <- Linear(DATAdum_ce)

### Intersect kernel PCA
pca_int <- kPCA(Kint,center=FALSE,plot=1:2,colors = "grey99")


#### Contribucions de cada al·lel als PCs
pcs_int <- kPCA_imp(DATAdum_ce,center = FALSE, secure=FALSE)

save(pcs_int$loadings,file="Results/PI_PCS_int.dat")
##SPECTRUM

load("Kernel_matrices/PI_mixSP1.dat")
load("Kernel_matrices/PI_mixSP2.dat")
load("Kernel_matrices/PI_mixSP3.dat")

## Spectrum kernel
### abund kmers
sp1_abund <- colMeans(mixSP1$feat_space)
sp2_abund <- colMeans(mixSP2$feat_space)
sp3_abund <- colMeans(mixSP3$feat_space)

mixSP2$feat_space <- desparsify(mixSP2$feat_space,2)
mixSP3$feat_space <- desparsify(mixSP3$feat_space,2)
dim(mixSP2$feat_space)
# 2212  346
dim(mixSP3$feat_space)
# 2212 1797

mixSP1$feat_space_ce <- scale(mixSP1$feat_space,center=T,scale=F)
mixSP2$feat_space_ce <- scale(mixSP2$feat_space,center=T,scale=F)
mixSP3$feat_space_ce <- scale(mixSP3$feat_space,center=T,scale=F)

### Càlcul kernel centrat
Ksp1 <- Linear(mixSP1$feat_space_ce)
Ksp2 <- Linear(mixSP2$feat_space_ce)
Ksp3 <- Linear(mixSP3$feat_space_ce)

### String kernel PCA
pca_s1 <- kPCA(Ksp1,center=F,plot=1:2,colors = "grey99")
pca_s2 <- kPCA(Ksp2,center=F,plot=1:2,colors = "grey99" )
pca_s3 <- kPCA(Ksp3,center=F,plot=1:2,colors = "grey99" )
# pca_s4 <- kPCA(Ksp4,center=F,plot=1:2,colors = "grey99" )

#### Contribucions de cada kmer als PCs.
pcs_s1 <- kPCA_imp(mixSP1$feat_space_ce,center = FALSE,projected=pca_s1$projection,secure=FALSE)
pcs_s2 <- kPCA_imp(mixSP2$feat_space_ce,center = FALSE,projected=pca_s2$projection,secure=FALSE)
pcs_s3 <- kPCA_imp(mixSP3$feat_space_ce,center = FALSE,projected=pca_s3$projection,secure=FALSE)


save(pcs_s1 ,file="Results/PI_PCS_mixSP1.dat")
save(pcs_s2 ,file="Results/PI_PCS_mixSP2.dat")
save(pcs_s3 ,file="Results/PI_PCS_mixSP3.dat")
