####################################
##### Mixtures vs no-mixtures ######
####################################

# Llibreries #
library(kerntools)
library(ggplot2)
library(cowplot)

# Declaració de funcions i constants #
source("helper_functions.R")


# Càrrega dels indexes #
load("Partitions/PI_IDX")
N <- length(c(IDX$work,IDX$additional_test_subtype))


# wild type B protease #
WTSeqseq <-toSeq(WTSeq)


# Càrrega de dades #

## Sense mescles d'aa
load("Kernel_matrices/PI_SP1_mean.dat")
load( ("Kernel_matrices/PI_SP2_mean.dat"))
load( ("Kernel_matrices/PI_SP3_mean.dat"))
load( ("Kernel_matrices/PI_Kov_mean.dat"))
load( ("Kernel_matrices/PI_LONG_15a_mean.dat"))
load("Kernel_matrices/PI_Kint.dat")
load("Kernel_matrices/PI_MIXLONG_15a.dat")
load("Kernel_matrices/PI_mixSP1.dat")
load("Kernel_matrices/PI_mixSP2.dat")
load("Kernel_matrices/PI_mixSP3.dat")


# Estudi previ propietats de les matrius de kernel #

### Distribució de valors
histK(cosNorm(Ksp1m),vn=FALSE, main = "PI Spectral l=1")
histK(cosNorm(Ksp2m),vn=FALSE, main = "PI Spectral l=2")
histK(cosNorm(mixSP2$K),vn=FALSE, main = "PI Spectral l=2 - MIXTURES")
histK(cosNorm(Ksp3m),vn=FALSE, main = "PI Spectral l=3")
histK(cosNorm(mixSP3$K),vn=FALSE, main = "PI Spectral l=3 - MIXTURES")

histK(cosNorm(Kovm),vn=FALSE, main = "PI Dirac")
histK(Kint$K,vn=FALSE, main = "PI Intersect")
histK(cosNorm(Long1m),vn=FALSE, main = "PI Interactions")
histK(cosNorm(Mixlong1$K),vn=FALSE, main = "PI Interactions-Intersect")


### Heatmaps, etc.  
heatK( (Ksp3m),raster = T,color = c("red","white"),cos.norm = T)
heatK( (mixSP3$K),raster = T,color = c("red","white"),cos.norm = T)
heatK( (Kovm),raster = T,color = c("red","white"),cos.norm = T)
heatK( (Kint$K),raster = T,color = c("red","white"),cos.norm = T)



# Kernel PCA results #

## Comparació dels kernels (3 string kernel i 1 Dirac kernel)
### Comparació de les PCA: coinèrcia
Ksim <- simK(list(Dirac=centerK(Kovm), Intersect = centerK(99*Kint$K),
                  Spectrum1=centerK(Ksp1m),mixSpectrum1=centerK(mixSP1$K),
                  Spectrum2=centerK(Ksp2m), mixSpectrum2=centerK(mixSP2$K),
                  Spectrum3=centerK(Ksp3m), mixSpectrum3=centerK(mixSP3$K),
                  Long1=centerK(Long1m), Mixlong1=centerK(Mixlong1$K)))
round(Ksim[seq(2,10,2),seq(2,10,2)],digits=3)


kPCA(Ksim,plot=1:2,y=rep(c("NoMixt","Mix"),5), labels=rep(c("Seq","1mers","2mers","3mers","long"),each=2))
diag(round(Ksim[seq(2,10,2),seq(1,10,2)],digits=3))


# ### Comparació de PCAs: Procrustes m
proc_pca <- matrix(0,nrow=1,ncol=5)
colnames(proc_pca) <- c("sp1","sp2","sp3","seq","long")
proc_pca["seq"] <- Procrustes(kPCA(Kovm,center=TRUE),kPCA(99*Kint$K,center=TRUE))$pro.cor
proc_pca["sp1"] <- Procrustes(kPCA(Ksp1m,center=TRUE),kPCA(mixSP1$K,center=TRUE))$pro.cor
proc_pca["sp2"] <- Procrustes(kPCA(Ksp2m,center=TRUE),kPCA(mixSP2$K,center=TRUE))$pro.cor
proc_pca["sp3"] <- Procrustes(kPCA(Ksp3m,center=TRUE),kPCA(mixSP3$K,center=TRUE))$pro.cor
proc_pca["long"] <- Procrustes(kPCA(Long1m,center=TRUE),kPCA(Mixlong1$K,center=TRUE))$pro.cor
proc_pca
 


Kovs <- list()

## Comparació amb els n=30 mostrejos
for(i in 1:30) {
  load(paste0("Kernel_matrices/PI_Kov_",i,".dat"))
  Kovs[[i]] <- centerK(Kov$K)
  rm(Kov)
  gc()
}
Kovs[[31]] <- centerK(Kovm)
Kovs[[32]] <- centerK(99*Kint$K)

Ksimov <- simK(Kovs)
Kintcenter <- kPCA(Ksimov,plot=1:2,y=factor(c(rep("",30),"mean","mixtures")))
Kintcenter$projection[31:32,] #la mitjana i intersect kernel estan quasi a l'origen de coordenades

## Comparació amb els n=30 mostrejos
Ksp3s <- list()
for(i in 1:30) {
  load(paste0("Kernel_matrices/PI_SP3_",i,".dat"))
  Ksp3s[[i]] <- centerK(SP3$K)
  rm(SP3)
  gc()
}
Ksp3s[[31]] <- centerK(Ksp3m)
Ksp3s[[32]] <- centerK(mixSP3$K)

Ksimov <- simK(Ksp3s)
summary(Ksimov[32,])
KSP3center <- kPCA(Ksimov,plot=1:2,y=factor(c(rep("",30),"mean","mixtures")))
KSP3center$projection[31:32,]



# PCS #


# Spectrum kernel
## length = 1

load("Results/PI_PCS_S1.dat")

nfeat <- 20
s1_firstpcs1 <- matrix("",nrow=30,ncol=nfeat)
for(i in 1:30) {
  s1_firstpcs1[i,] <-   plotImp(PCS_S1[[i]][1,], nfeat=nfeat,  ylegend="PC2",absolute=FALSE, relative = FALSE,
                                main="PC1", leftmargin = 5,  color="pink")$first_features
}
summary(toFactor(s1_firstpcs1))

## length = 2
load("Results/PI_PCS_S2.dat")
nfeat <- 20
s2_firstpcs1 <- matrix("",nrow=30,ncol=nfeat)
for(i in 1:30) {
  s2_firstpcs1[i,] <-   plotImp(PCS_S2[[i]][1,], nfeat=nfeat,  ylegend="PC2",absolute=FALSE, relative = FALSE,
                                main="PC1", leftmargin = 5,  color="pink")$first_features
}
summary(toFactor(s2_firstpcs1))
#els top20 amb mixtures estan tots a no-mixtures

## length = 3
load("Results/PI_PCS_S3.dat")
nfeat <- 20
s3_firstpcs1 <- s3_firstpcs2 <- matrix("",nrow=30,ncol=nfeat)
for(i in 1:30) {
  s3_firstpcs1[i,] <-   plotImp(PCS_S3[[i]][1,], nfeat=nfeat, absolute=FALSE, relative = FALSE,
                                main="PC1", leftmargin = 5,  color="pink")$first_features
  s3_firstpcs2[i,] <-   plotImp(PCS_S3[[i]][2,], nfeat=nfeat,  absolute=FALSE, relative = FALSE,
                                main="PC2", leftmargin = 5,  color="pink")$first_features
}
summary(toFactor(s3_firstpcs1))
summary(toFactor(s3_firstpcs2))

sp3_pc1 <- plotImp(pcs_s3$loadings[1,], y=pcs_s3$loadings[2,], ylegend="PC2",nfeat=20,absolute=FALSE, relative = FALSE,
                   main="PC1", leftmargin = 5,  color="pink")
setdiff(unique(as.vector(s3_firstpcs1)),sp3_pc1$first_features) #només falta MTQ (que està immediatament després)



## Amb mescles d'aa:

load(file="Results/PI_PCS_mixSP1.dat")
load(file="Results/PI_PCS_mixSP2.dat")
load(file="Results/PI_PCS_mixSP3.dat")

pcs_s1$loadings[1:2,]
pcs_s2$loadings[1:2,]
pcs_s3$loadings[1:2,]

load("Results/PI_PCS_S1.dat")
pcs_s1_comp <- matrix(0,nrow=31,ncol=length(AA))
colnames(pcs_s1_comp) <- AA
for(i in 1:30) pcs_s1_comp[i,colnames(PCS_S1[[i]])] <- PCS_S1[[i]][1,]
pcs_s1_comp[31,colnames(pcs_s1$loadings)] <- pcs_s1$loadings[1,] ## això ve s'altre document "no mixtures"
pcs_s1_comp <- desparsify(pcs_s1_comp,dim=1:2)
sort(colMeans(abs(pcs_s1_comp)),decreasing = TRUE)
Kken <- Kendall(pcs_s1_comp,samples.in.rows = TRUE)
heatK(Kken)
histK(Kken)
summary(Kken[31,-31])


load("Results/PI_PCS_S2.dat")
pcs_s2_comp <- matrix(0,nrow=31,ncol=length(AA)^2)
colnames(pcs_s2_comp) <- permute_rep(alphabet=AA,l=2)
for(i in 1:30) pcs_s2_comp[i,colnames(PCS_S2[[i]])] <- PCS_S2[[i]][1,]
pcs_s2_comp[31,colnames(pcs_s2$loadings)] <- pcs_s2$loadings[1,] ## això ve s'altre document "no mixtures"
pcs_s2_comp <- desparsify(pcs_s2_comp,dim=1:2)
sort(colMeans(abs(pcs_s2_comp)),decreasing = TRUE)
Kken <- Kendall(pcs_s2_comp,samples.in.rows = TRUE)
heatK(Kken)
histK(Kken,breaks=100)
summary(Kken[31,-31])
kPCA(K=Kken,plot=1:2,y = as.factor(c(rep("No Mixtures",30),"Mixtures")))


pcs_sp3_comp <- matrix(0,nrow=31,ncol=ncol(pcs_s3$loadings))
colnames(pcs_sp3_comp) <- colnames(pcs_s3$loadings)
for(i in 1:30) pcs_sp3_comp[i,colnames(PCS_S3[[i]])] <- PCS_S3[[i]][1,]
pcs_sp3_comp[31,colnames(pcs_s3$loadings)] <- pcs_s3$loadings[1,]
pcs_sp3_comp <- desparsify(pcs_sp3_comp,dim=1:2)
sort(colMeans(abs(pcs_sp3_comp)),decreasing = TRUE)
Kken3 <- Kendall(pcs_sp3_comp,samples.in.rows = TRUE)
Kken3
histK(Kken3,breaks=100)
summary(Kken3[31,-31])
panelB <- kPCA(K=Kken3,plot=1:2,title="Spectrum l=3",y = as.factor(c(rep("No Mixtures",30),"Mixtures")))
legend <- cowplot::get_legend(panelB$plot)
legend <- ggplotGrob(panelB$plot + theme(legend.position = "bottom"))$grobs
legend <- legend[[which(sapply(legend, function(x) x$name) == "guide-box")]]

# Dirac kernel
load("Results/PI_PCS_ov.dat")
load("Results/PI_PCS_int.dat")
nfeat <- 20
ov_firstpcs1 <- matrix("",nrow=30,ncol=nfeat)
for(i in 1:30) {
  ov_firstpcs1[i,] <-   plotImp(PCS_OV[[i]][1,], nfeat=nfeat,  ylegend="PC2",absolute=FALSE, relative = FALSE,
                                main="PC1", leftmargin = 5,  color="pink")$first_features
}
summary(toFactor(ov_firstpcs1))
setdiff(unique(as.vector(ov_firstpcs1)),
        int_firstpcs1 <-plotImp(pcs_int$loadings[1,], nfeat=nfeat,  ylegend="PC2",absolute=FALSE, relative = FALSE,
                                main="PC1", leftmargin = 5,  color="pink")$first_features)

pcs_ov_comp <- matrix(0,nrow=31,ncol=length(PrPositions)*20)
colnames(pcs_ov_comp) <- paste(rep(PrPositions,each=20),rep(AA,99),sep="_")
for(i in 1:30) pcs_ov_comp[i,colnames(PCS_OV[[i]])] <- PCS_OV[[i]][1,]
pcs_ov_comp[31,colnames(pcs_int$loadings)] <- pcs_int$loadings[1,] ## això ve s'altre document "no mixtures"
pcs_ov_comp <- desparsify(pcs_ov_comp,dim=1:2)
sort(colMeans(abs(pcs_ov_comp)),decreasing = TRUE)
Kken <- Kendall(pcs_ov_comp,samples.in.rows = TRUE)
heatK(Kken)
summary(Kken[31,-31])
panelA <- kPCA(K=Kken,plot=1:2,title="Dirac vs Intersect",y =  as.factor(c(rep("No Mixtures",30),"Mixtures")))$plot + ggplot2::theme(legend.position="none")

plot_grid(plot_grid(panelA,panelB$plot + ggplot2::theme(legend.position="none"),
                    labels = c("A","B"),ncol=2),
          legend, nrow=2, rel_heights = c(0.9,0.1))
