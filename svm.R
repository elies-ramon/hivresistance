#################
##### SVMs ######
#################


# Llibreries #
library(kerntools)
library(kernlab)

# Declaració de constants i funcions #
source("helper_functions.R")

# Càrrega dels indexes #
load("Partitions/PI_IDX")


# Càrrega de les dades. Particionat #

load("Kernel_matrices/PI_Kint.dat")
load("Kernel_matrices/PI_mixSP1.dat")
load("Kernel_matrices/PI_mixSP2.dat")
load("Kernel_matrices/PI_mixSP3.dat")
load("Kernel_matrices/PI_MIXLONG1.dat")

### llevar interaccions dins d'una mateixa pos
for(pos in PrPositions) {
  delete <- paste(rep(pos,2),collapse ="_")
  Mixlong1$feat_space[,grep(delete,colnames(Mixlong1$feat_space))] <- 0
}

Mixlong1$feat_space <- desparsify(Mixlong1$feat_space)
dim(Mixlong1$feat_space)
head(colnames(Mixlong1$feat_space))

Mixlong1$K <- Linear(Mixlong1$feat_space)
histK( cosNorm(Mixlong1$K ))

load("Data/PI_MultiA_clean")
rownames(DATA) <- 1:nrow(DATA)

te_add_subtype <-  DATA[IDX$additional_test_subtype,"SeqID"]
te_idx <- DATA[ IDX$test,"SeqID"]


DATA <- droplevels(DATA[as.numeric(c(IDX$work,IDX$additional_test_subtype)),])
dim(DATA)  #2212  115

Partition <- rep("Training",nrow(DATA))
Partition[which( DATA$SeqID %in% te_idx)] <- "Test"
Partition[which( DATA$SeqID %in% te_add_subtype)] <- "Additional test"
summary(factor(Partition))
# Additional test            Test        Training
#            125             407            1680


## Comprovacions, idx
summary(factor(DATA$Method))
summary(factor(DATA$Type[Partition!="Additional test"]))
sum(DATA[Partition=="Additional test","Subtype"] == "B")

anyDuplicated(te_idx)
anyDuplicated(te_add_subtype)

te_idx <-  which(Partition == "Test")
train_idx <-  which(Partition == "Training")
alt_te_idx <-  which(Partition == "Additional test")

intersect(DATA[ train_idx,"PtID"],DATA[ te_idx,"PtID"])
intersect(te_idx,train_idx)
length(te_idx) + length(train_idx) + length(IDX$additional_test_subtype) == nrow(DATA)



# Entrenament #

nmodels <- 5

## 5x2 CV

  val_idx <- list()
  val_N <- round(length(train_idx)/5)
  sampled_train <- sample(train_idx)
  for(i in 0:4)   val_idx[[i+1]] <- sampled_train[(1+(i*val_N)):((1+i)*val_N)]
  val_idx[[i+1]] <- val_idx[[i+1]] [!is.na( val_idx[[i+1]] )]


  C <- c(0.1,1,10,100) # hiperparàmetre C
  HYP <- array(0,dim=c(length(C),ncol=length(DRUGS),nmodels))
  rownames(HYP) <- C

  FIT <-  array(NA,dim=c(nrow(DATA),ncol=length(DRUGS),nmodels))
  colnames(HYP) <- colnames(FIT) <- DRUGS
  VALS <- PREDS <- FIT

  SV <- array(0,dim=c(nrow(DATA),length(DRUGS),nmodels),
              dimnames = list(1:nrow(DATA),DRUGS,c("Spectrum1","Spectrum2","Spectrum3","Intersect","Interactions")))


  ### Importances
  #### Intersect
  #Reorganitzar + Esborrar posicions que apareixen 0 voltes:
  DATAdum_delete <- rearrange_kint(Kint$feat_space,rows=rownames(DATA),AA=AA,PrPositions = PrPositions)

  Kint_imps <- matrix(0,nrow= ncol(DATAdum_delete),ncol=length(DRUGS))
  colnames(Kint_imps) <- DRUGS
  rownames(Kint_imps) <- colnames(DATAdum_delete)

  ####SPECTRUM L=1
  Sp1_imps <- matrix(0,nrow= ncol(mixSP1$feat_space),ncol=length(DRUGS))
  colnames(Sp1_imps) <- DRUGS
  rownames(Sp1_imps) <- colnames(mixSP1$feat_space)

  ####SPECTRUM L=2
  Sp2_imps <- matrix(0,nrow= ncol(mixSP2$feat_space),ncol=length(DRUGS))
  colnames(Sp2_imps) <- DRUGS
  rownames(Sp2_imps) <- colnames(mixSP2$feat_space)

  ####SPECTRUM L=3
  Sp3_imps <- matrix(0,nrow= ncol(mixSP3$feat_space),ncol=length(DRUGS))
  colnames(Sp3_imps) <- DRUGS
  rownames(Sp3_imps) <- colnames(mixSP3$feat_space)

  ####mixlong
  DATAdum_mix <- desparsify(Mixlong1$feat_space)

  Mixlong1_imps <- matrix(0,nrow= ncol(DATAdum_mix),ncol=length(DRUGS))
  colnames(Mixlong1_imps) <- DRUGS
  rownames(Mixlong1_imps) <- colnames(DATAdum_mix)

  ## Normalització cosinus (valors entre 0 i 1)
  Kint <- cosNorm(Kint$K)
  Ksp1 <- cosNorm(mixSP1$K)
  Ksp2 <- cosNorm(mixSP2$K)
  Ksp3 <- cosNorm(mixSP3$K)
  Mixlong1 <- cosNorm(Mixlong1$K)

for(i in 1:5) {
  for(drug in DRUGS) {
    print(drug)
    Y.tr <- DATA[ ,drug]
    train_drug <- setdiff(train_idx, which(is.na(Y.tr)))
    train_drug <- setdiff(train_drug, val_idx[[i]])
    Y.tr <- Y.tr[train_drug]
    Kov.train <- Kint[train_drug,train_drug]
    Ksp1.train <- Ksp1[train_drug,train_drug]
    Ksp2.train <- Ksp2[train_drug,train_drug]
    Ksp3.train <- Ksp3[train_drug,train_drug]
    Long1.train <- Mixlong1[train_drug,train_drug]


    # Grid search 2-cross-validation
    for(c in C) {
      ksp1_model <- kernlab::ksvm(Ksp1.train,Y.tr,type="eps-svr",kernel="matrix",cross=2,C=c)
      HYP[as.character(c),drug,1] <- kernlab::cross(ksp1_model)

      ksp2_model <- kernlab::ksvm(Ksp2.train,Y.tr,type="eps-svr",kernel="matrix",cross=2,C=c)
      HYP[as.character(c),drug,2] <- kernlab::cross(ksp2_model)

      ksp3_model <- kernlab::ksvm(Ksp3.train,Y.tr,type="eps-svr",kernel="matrix",cross=2,C=c)
      HYP[as.character(c),drug,3] <- kernlab::cross(ksp3_model)

      kov_model <- kernlab::ksvm(Kov.train,Y.tr,type="eps-svr",kernel="matrix",cross=2,C=c)
      HYP[as.character(c),drug,4] <- kernlab::cross(kov_model)

      long1_model <- kernlab::ksvm(Long1.train,Y.tr,type="eps-svr",kernel="matrix",cross=2,C=c)
      HYP[as.character(c),drug,5] <- kernlab::cross(long1_model)

    }

    # Millor model
    ksp1_best <- kernlab::ksvm(Ksp1.train,Y.tr, scaled=FALSE,
                               type="eps-svr",kernel="matrix",
                               C=as.numeric(names(which.min(HYP[,drug,1])))) # Rular el mètode

    ksp2_best <- kernlab::ksvm(Ksp2.train,Y.tr, scaled=FALSE,
                               type="eps-svr",kernel="matrix",
                               C=as.numeric(names(which.min(HYP[,drug,2])))) # Rular el mètode

    ksp3_best <- kernlab::ksvm(Ksp3.train,Y.tr, scaled=FALSE,
                               type="eps-svr",kernel="matrix",
                               C=as.numeric(names(which.min(HYP[,drug,3])))) # Rular el mètode

    kov_best <- kernlab::ksvm(Kov.train,Y.tr, scaled=FALSE,
                              type="eps-svr",kernel="matrix",
                              C=as.numeric(names(which.min(HYP[,drug,4]))))

    long1_best <- kernlab::ksvm(Long1.train,Y.tr, scaled=FALSE,
                                type="eps-svr",kernel="matrix",
                                C=as.numeric(names(which.min(HYP[,drug,5])))) # Rular el mètode


    ## Calcula la matriu de test, que en realitat és una matriu test vs SV
    ### Primer s'ha de calcular la matriu test vs train
    #
    Kov.te <- Kint[val_idx[[i]],train_drug]
    Ksp1.te <- Ksp1[val_idx[[i]],train_drug]
    Ksp2.te <- Ksp2[val_idx[[i]],train_drug]
    Ksp3.te <- Ksp3[val_idx[[i]],train_drug]
    Long1.te <- Mixlong1[val_idx[[i]],train_drug]

    Kov.te <- Kov.te[,SVindex(kov_best),drop=FALSE]
    Kov.te <- as.kernelMatrix(Kov.te)

    Ksp1.te <- Ksp1.te[,SVindex(ksp1_best),drop=FALSE]
    Ksp1.te <- as.kernelMatrix(Ksp1.te)

    Ksp2.te <- Ksp2.te[,SVindex(ksp2_best),drop=FALSE]
    Ksp2.te <- as.kernelMatrix(Ksp2.te)

    Ksp3.te <- Ksp3.te[,SVindex(ksp3_best),drop=FALSE]
    Ksp3.te <- as.kernelMatrix(Ksp3.te)

    Long1.te <- Long1.te[,SVindex(long1_best),drop=FALSE]
    Long1.te <- as.kernelMatrix(Long1.te)



    VALS[val_idx[[i]],drug,4] <- kernlab::predict(kov_best,Kov.te) ## Valors predits
    VALS[val_idx[[i]],drug,1] <- kernlab::predict(ksp1_best,Ksp1.te) ## Valors predits
    VALS[val_idx[[i]],drug,2] <- kernlab::predict(ksp2_best,Ksp2.te) ## Valors predits
    VALS[val_idx[[i]],drug,3] <- kernlab::predict(ksp3_best,Ksp3.te) ## Valors predits
    VALS[val_idx[[i]],drug,5] <- kernlab::predict(long1_best,Long1.te) ## Valors predits

  }
}

  #Guardam dades

  save(val_idx,file="Results/val_idx.dat")
  save(VALS,file="Results/validation_results.dat")
  save(HYP,file="Results/HYP_results.dat")



 ## Models definitius

for(drug in DRUGS) {
  print(drug)
  Y.tr <- DATA[ ,drug]
  train_drug <- setdiff(train_idx, which(is.na(Y.tr)))
  te_drug <- c(te_idx, alt_te_idx)
  Y.tr <- Y.tr[train_drug]
  Kint.train <- Kint[train_drug,train_drug]
  Ksp1.train <- Ksp1[train_drug,train_drug]
  Ksp2.train <- Ksp2[train_drug,train_drug]
  Ksp3.train <- Ksp3[train_drug,train_drug]
  Mixlong1.train <- Mixlong1[train_drug,train_drug]

  # Millor model
  ksp1_best <- kernlab::ksvm(Ksp1.train,Y.tr, scaled=FALSE,
                             type="eps-svr",kernel="matrix",
                             C=as.numeric(names(which.min(HYP[,drug,1])))) # Rular el mètode

  ksp2_best <- kernlab::ksvm(Ksp2.train,Y.tr, scaled=FALSE,
                             type="eps-svr",kernel="matrix",
                             C=as.numeric(names(which.min(HYP[,drug,2])))) # Rular el mètode

  ksp3_best <- kernlab::ksvm(Ksp3.train,Y.tr, scaled=FALSE,
                             type="eps-svr",kernel="matrix",
                             C=as.numeric(names(which.min(HYP[,drug,3])))) # Rular el mètode

  kint_best <- kernlab::ksvm(Kint.train,Y.tr, scaled=FALSE,
                            type="eps-svr",kernel="matrix",
                            C=as.numeric(names(which.min(HYP[,drug,4]))))

  mixlong1_best <- kernlab::ksvm(Mixlong1.train,Y.tr, scaled=FALSE,
                              type="eps-svr",kernel="matrix",
                              C=as.numeric(names(which.min(HYP[,drug,5])))) # Rular el mètode


  #Guardam fitted values
  FIT[train_drug,drug,1] <- fitted(kint_best)
  FIT[train_drug,drug,2] <- fitted(ksp1_best)
  FIT[train_drug,drug,3] <- fitted(ksp2_best)
  FIT[train_drug,drug,4] <- fitted(ksp3_best)
  FIT[train_drug,drug,5] <- fitted(mixlong1_best)


  #Guardam vectors suport
  SV[train_drug[alphaindex(ksp1_best)],drug,1] <- 1
  SV[train_drug[alphaindex(ksp2_best)],drug,2] <- 1
  SV[train_drug[alphaindex(ksp3_best)],drug,3] <- 1
  SV[train_drug[alphaindex(kint_best)],drug,4] <- 1
  SV[train_drug[alphaindex(mixlong1_best)],drug,5] <- 1


  ## Calcula la matriu de test, que en realitat és una matriu test vs SV
  ### Primer s'ha de calcular la matriu test vs train

  Kint.te <- Kint[te_idx,train_drug]
  Ksp1.te <- Ksp1[te_idx,train_drug]
  Ksp2.te <- Ksp2[te_idx,train_drug]
  Ksp3.te <- Ksp3[te_idx,train_drug]
  Mixlong1.te <- Mixlong1[te_idx,train_drug]

  Kint.te <- Kint.te[,SVindex(kint_best),drop=FALSE]
  Kint.te <- as.kernelMatrix(Kint.te)

  Ksp1.te <- Ksp1.te[,SVindex(ksp1_best),drop=FALSE]
  Ksp1.te <- as.kernelMatrix(Ksp1.te)

  Ksp2.te <- Ksp2.te[,SVindex(ksp2_best),drop=FALSE]
  Ksp2.te <- as.kernelMatrix(Ksp2.te)

  Ksp3.te <- Ksp3.te[,SVindex(ksp3_best),drop=FALSE]
  Ksp3.te <- as.kernelMatrix(Ksp3.te)

  Mixlong1.te <- Mixlong1.te[,SVindex(mixlong1_best),drop=FALSE]
  Mixlong1.te <- as.kernelMatrix(Mixlong1.te)

  PREDS[te_idx,drug,1] <- kernlab::predict(ksp1_best,Ksp1.te) ## Valors predits
  PREDS[te_idx,drug,2] <- kernlab::predict(ksp2_best,Ksp2.te) ## Valors predits
  PREDS[te_idx,drug,3] <- kernlab::predict(ksp3_best,Ksp3.te) ## Valors predits
  PREDS[te_idx,drug,4] <- kernlab::predict(kint_best,Kint.te) ## Valors predits
  PREDS[te_idx,drug,5] <- kernlab::predict(mixlong1_best,Mixlong1.te) ## Valors predits

  #Test alternatiu
  Kint.tealt <- Kint[alt_te_idx,train_drug]
  Ksp1.tealt <- Ksp1[alt_te_idx,train_drug]
  Ksp2.tealt <- Ksp2[alt_te_idx,train_drug]
  Ksp3.tealt <- Ksp3[alt_te_idx,train_drug]
  Mixlong1.tealt <- Mixlong1[alt_te_idx,train_drug]

  Kint.tealt <- Kint.tealt[,SVindex(kint_best),drop=FALSE]
  Kint.tealt <- as.kernelMatrix(Kint.tealt)

  Ksp1.tealt <- Ksp1.tealt[,SVindex(ksp1_best),drop=FALSE]
  Ksp1.tealt <- as.kernelMatrix(Ksp1.tealt)

  Ksp2.tealt <- Ksp2.tealt[,SVindex(ksp2_best),drop=FALSE]
  Ksp2.tealt <- as.kernelMatrix(Ksp2.tealt)

  Ksp3.tealt <- Ksp3.tealt[,SVindex(ksp3_best),drop=FALSE]
  Ksp3.tealt <- as.kernelMatrix(Ksp3.tealt)

  Mixlong1.tealt <- Mixlong1.tealt[,SVindex(mixlong1_best),drop=FALSE]
  Mixlong1.tealt <- as.kernelMatrix(Mixlong1.tealt)


  PREDS[alt_te_idx,drug,1] <- kernlab::predict(ksp1_best,Ksp1.tealt) ## Valors predits
  PREDS[alt_te_idx,drug,2] <- kernlab::predict(ksp2_best,Ksp2.tealt) ## Valors predits
  PREDS[alt_te_idx,drug,3] <- kernlab::predict(ksp3_best,Ksp3.tealt) ## Valors predits
  PREDS[alt_te_idx,drug,4] <- kernlab::predict(kint_best,Kint.tealt) ## Valors predits
  PREDS[alt_te_idx,drug,5] <- kernlab::predict(mixlong1_best,Mixlong1.tealt) ## Valors predits


  # Guardam importàncies
  Imps <- svm_imp(DATAdum_delete[train_drug,],alphaindex(kint_best) ,coef(kint_best),result = "original" ,scale = F,center = F,cos.norm = T)
  Kint_imps[names(Imps),drug] <- Imps

  Imps <- svm_imp(mixSP1$feat_space[train_drug,],alphaindex(ksp1_best) ,coef(ksp1_best),result = "original" ,scale = F,center = F,cos.norm = T)
  Sp1_imps[names(Imps),drug] <- Imps

  Imps <- svm_imp(mixSP2$feat_space[train_drug,],alphaindex(ksp2_best) ,coef(ksp2_best),result = "original" ,scale = F,center = F,cos.norm = T)
  Sp2_imps[names(Imps),drug] <- Imps

  Imps <- svm_imp(mixSP3$feat_space[train_drug,],alphaindex(ksp3_best) ,coef(ksp3_best),result = "original" ,scale = F,center = F,cos.norm = T)
  Sp3_imps[names(Imps),drug] <- Imps

  Imps <- svm_imp(DATAdum_mix[train_drug,],alphaindex(mixlong1_best) ,coef(mixlong1_best),result = "original" ,scale = F,center = F,cos.norm = T)
  Mixlong1_imps[names(Imps),drug] <- Imps
}



# Guardar resultats

IMPS <- list(Spectrum1=Sp1_imps,Spectrum2=Sp2_imps,Spectrum3=Sp3_imps,Intersect=Kint_imps,Interactions=Mixlong1_imps)
save(PREDS,file="Results/PI_PREDS.dat")
save(FIT,file="Results/PI_FIT.dat")
save(IMPS, file="Results/PI_IMPS.dat")
save(SV, file="Results/PI_SV.dat")

