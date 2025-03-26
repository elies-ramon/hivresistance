##############################################
##### PAIRWISE DISTANCES BETWEEN RESIDUES#####
##############################################

# Llibreries #
library(kerntools)
library(ggplot2)

source("helper_functions.R")

# Càrrega de dades #

#### PDB 3OXC per a HIV-1 wt proteasa (com està descrit a Xiaxia Yu, 2014)
structure <- read.table("Data/3oxc.proc2",header = TRUE,stringsAsFactors = FALSE)
p <- cbind(structure$X,structure$Y,structure$Z)
colnames(p) <- c("X","Y","Z")
rownames(p) <- paste("P",structure$Pos.in.subunit,sep="")
d <- nrow(p)


# Distància #
## La proteasa és un homodimer. Distancia entre carbonis alfa:
dist.ind <- expand.grid.mod(1:d,rep=FALSE)
DIST <- matrix(0,nrow=d,ncol=d)
for(i in 1:nrow(dist.ind)) {
  a <- dist.ind[i,1]
  b <- dist.ind[i,2]
  DIST[a,b] <- DIST[b,a] <- calc_dist(p[a,],p[b,])
}

## Prenem la mínima distància possible:
DIST99  <- list()
DIST99[[1]] <- DIST[1:99,1:99]
DIST99[[2]]  <- DIST[1:99,100:198]
DIST99[[3]] <-  DIST[100:198,100:198]
DIST99 <-  Reduce(pmin, DIST99)

colnames(DIST99) <- rownames(DIST99) <- PrPositions
summary(as.vector(DIST99))

save(DIST99,file="LD/PI_DIST.dat")

heatK(DIST99,raster = F,title = "Minimum distance between residues (in angstrom)") + scale_fill_gradientn(colours=c("firebrick1","white","turquoise1"))


