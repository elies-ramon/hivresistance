##################################
##### LINKAGE DISEQUILIBRIUM ######
####################################

# Llibreries #
library(kerntools)
library(ggplot2)

source("helper_functions.R")


# Càrrega de dades #


# Càrrega dels indexes #
load("Partitions/PI_IDX")

load("Data/PI_Seq_clean")
LengthP <- length(PrPositions)
LD <- array(0,dim = c(LengthP,LengthP,30))
colnames(LD) <- rownames(LD) <- PrPositions

Rsquared <-list()

for (i in 1:30)   {
  DATA.l <-   DATA.list[[i]]

  # Particionat
  DATA <- droplevels(DATA.l[c(IDX$work,IDX$additional_test_subtype),])
  dim(DATA)  #2212  115

  ### LD
  #FREQ0 = cada entrada són les freq al·lèliques d'una posició de la proteïna.
  FREQ0 <- calc_freq(DATA[,PrPositions])
  FREQ0 <- unlist(FREQ0)
  names(FREQ0) <- sub("\\.","_",names(FREQ0))

  #FREQ1 = cada entrada és P1P2...Pd-1Pd, les mult de freq alel de tots els al·lels de les 2 proteïnes
  FREQ1 <- combn(FREQ0,2,FUN=prod)
  cols1 <- t(combn(sub(".*_", "", names(FREQ0)),2))
  cols2 <- t(combn(sub("_.*", "", names(FREQ0)),2))
  names(FREQ1) <- paste(paste(cols2[,1],cols2[,2],sep="_"), paste0(cols1[,1],cols1[,2]) ,sep="_")
  ### més variants necessàries per a calcular Dmax
  FREQ1minus <- combn((1-FREQ0),2,FUN=prod)
  # FREQ1minus1 <- combn(FREQ0,2,FUN=function(x)return(x[1]*(1-x[2])))
  # FREQ1minus2 <- combn(FREQ0,2,FUN=function(x)return((1-x[1])*x[2]))
  # names(FREQ1minus) <-  names(FREQ1minus1) <-  names(FREQ1minus2) <-  names(FREQ1)
  names(FREQ1minus) <-  names(FREQ1)

  ## retirar els que venen de sa mateixa posició
  deletethis <- paste0(paste(matrix(PrPositions),matrix(PrPositions),sep="_"),"_")
  for(k in deletethis)  FREQ1 <- FREQ1[ !grepl(k,names(FREQ1))]
  FREQ1minus <- FREQ1minus[  names(FREQ1)]
  # FREQ1minus1 <- FREQ1minus1[names(FREQ1)]
  # FREQ1minus2 <- FREQ1minus2[ names(FREQ1)]
  length(FREQ1minus)

  #FREQ2 = cada entrada és P12...Pd-1d, les freq alel observades en la realitat
  DATAdum2 <- apply(apply(DATA[,PrPositions],2,as.character),1,pos_interac)
  DATAdum2 <- t(DATAdum2)
  cols <- t(combn(PrPositions,2))
  colnames(DATAdum2) <- paste(cols[,1],cols[,2],sep="_")
  FREQ2 <-  calc_freq(toFactor(DATAdum2))
  FREQ2 <- unlist(FREQ2)
  names(FREQ2) <- sub("\\.","_",names(FREQ2))
  length(FREQ2)
  cols3 <- names(FREQ1)[!(names(FREQ1) %in% names(FREQ2))]
  notin <- rep(0,length(cols3))
  names(notin) <- cols3
  FREQ2  <-  c(FREQ2 , notin)
  length(FREQ2)

  FREQ2 <- FREQ2[names(FREQ1)]
  D <- FREQ2-FREQ1
  rsquared <- D^2 / (FREQ1minus)
  summary(rsquared)
  rsquared[is.nan(rsquared)] <- 0
  Rsquared[[i]] <- rsquared/FREQ1

  sort(rsquared,decreasing = T)[1:100]

  for(k in 1:(LengthP-1)) {
    for(j in (k+1):LengthP) {
      A <- PrPositions[k]
      B <- PrPositions[j]
     LD[A,B,i] <- sum(rsquared[grep(paste0(A, "_",B,"_"),names(rsquared))])
    }
  }
  print(heatK(t(LD[,,i]),cos.norm = FALSE,color = c("white","red"),title=paste("Protease",i,"LD between protein positions")))
  LD[,,i] <- LD[,,i] + t(LD[,,i])
  diag(LD[,,i]) <- 1
}
save(LD,file="LD/PI_LD.dat")
save(Rsquared,file="LD/PI_Rsqd.dat")

print(heatK(rowMeans(LD,dims = 2),cos.norm = FALSE,color = c("white","red"),title=paste("Protease average LD between protein positions")))
sort(LD,decreasing = T)


######### AMB LES MESCLES D'AA

# Càrrega dels indexes #
load("Partitions/PI_IDX")

# Càrrega de dades #
load("Data/PI_MultiA_clean")
rownames(DATA) <- 1:nrow(DATA)

DATA <- droplevels(DATA[as.numeric(c(IDX$work,IDX$additional_test_subtype)),])
dim(DATA)  #2212  115


Kint <- Intersect (DATA[,PrPositions],elements = AA,comp = "sum" ,feat_space=TRUE )
FREQ0_mixt <- rearrange_kint(Kint$feat_space,rows=rownames(DATA),AA=AA,PrPositions = PrPositions)
mixtures <- aggregate_imp(FREQ0_mixt)
rownames(mixtures) <- sub("X","",rownames(mixtures))
mixtures_where <- which(mixtures>1,arr.ind = T)
for(j in unique(mixtures_where[,2])) {
  samples <- mixtures_where[mixtures_where[,2]==j,1]
  pos <- grep(paste0("^P",j,"_"),colnames(FREQ0_mixt),value=TRUE)
  FREQ0_mixt[samples,pos] <- FREQ0_mixt[samples,pos]/mixtures[samples,j]
 }
 sort(colMeans(FREQ0_mixt)) ## freq amb els al·lels en mescla.
 save(FREQ0_mixt,file="LD/Freqs_mixt.dat")

 LengthP <- length(PrPositions)
 LD <- array(0,dim = c(LengthP,LengthP))
 colnames(LD) <- rownames(LD) <- PrPositions

   ### LD
   #FREQ0 = cada entrada són les freq al·lèliques d'una posició de la proteïna.
  FREQ0 <- colMeans(FREQ0_mixt)

   #FREQ1 = cada entrada és P1P2...Pd-1Pd, les mult de freq alel de tots els al·lels de les 2 proteïnes
   FREQ1 <- combn(FREQ0,2,FUN=prod)
   cols1 <- t(combn(sub(".*_", "", names(FREQ0)),2))
   cols2 <- t(combn(sub("_.*", "", names(FREQ0)),2))
   names(FREQ1) <- paste(paste(cols2[,1],cols2[,2],sep="_"), paste0(cols1[,1],cols1[,2]) ,sep="_")
   ### més variants necessàries per a calcular Dmax
   FREQ1minus <- combn((1-FREQ0),2,FUN=prod)
   names(FREQ1minus) <-  names(FREQ1)

   ## retirar els que venen de sa mateixa posició
   deletethis <- paste0(paste(matrix(PrPositions),matrix(PrPositions),sep="_"),"_")
   for(k in deletethis)  FREQ1 <- FREQ1[ !grepl(k,names(FREQ1))]
   FREQ1minus <- FREQ1minus[  names(FREQ1)]
   length(FREQ1minus)

   #FREQ2 = cada entrada és P12...Pd-1d, les freq alel observades en la realitat
   DATAdum_delete <- desparsify(FREQ0_mixt,2)
   dim(DATAdum_delete)
   DATAlong1 <- matrix(0,nrow=nrow(DATAdum_delete),ncol=choose(ncol(DATAdum_delete),2))
   for(i in 1:nrow(DATAdum_delete)) {
     DATAlong1[i,] <- combn(DATAdum_delete[i,],2,FUN=prod)
   }
   cols <- t(combn(colnames(DATAdum_delete),2))
   cols1 <-  sub(".*_", "", cols)
   cols2 <-  sub("_.*", "", cols)
   colnames(DATAlong1) <- paste(paste(cols2[,1],cols2[,2],sep="_"), paste0(cols1[,1],cols1[,2]) ,sep="_")
   FREQ2 <-  colMeans(DATAlong1)
   length(FREQ2)

   FREQ2 <- FREQ2[names(FREQ1)]
   D <- FREQ2-FREQ1

   rsquared <- D^2 / (FREQ1minus[names(D)])
   summary(rsquared)
    rsquared[is.nan(rsquared)] <- 0
   Rsquared  <-  rsquared/FREQ1

   sort(Rsquared,decreasing = T)[1:100]

   for(k in 1:(LengthP-1)) {
     for(j in (k+1):LengthP) {
       A <- PrPositions[k]
       B <- PrPositions[j]
       LD[A,B] <- sum(rsquared[grep(paste0(A, "_",B,"_"),names(rsquared))])
     }
   }
   max(LD)
   LD  <- LD + t(LD )
   heatK( LD,cos.norm = FALSE,color = c("white","red"),title= "Protease LD between protein positions")

 save(LD,file="LD/PI_LD.dat")
 save(Rsquared,file="LD/PI_Rsqd.dat")

 ## Estudi LD
 LD <- rowMeans(LD,dims = 2)
 diag(LD) <- NA

 sort(Rsquared[[1]],decreasing = T)[1:20]


 sort(LD,decreasing = T)[1:20]
 first10 <- unique(sort(LD,decreasing = T)[1:20])
 for(pos in first10) {
   print(which(LD==pos,arr.ind = T))
 }
 # 25-32, 47-53, 76 and 80-84
  
 heatK(as.matrix(LD),title = "Linkage Disequilibrium",raster = F,color = c("grey99","firebrick2")) +
   theme(  axis.text.x = element_text(size = 8),  axis.text.y= element_text(size = 8) )

 mostassoc <-  sort(names(sort(rowMeans(LD,na.rm = T),decreasing = T)[1:10]))
 mostassoc
 