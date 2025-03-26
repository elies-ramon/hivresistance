###############################
##### FUNCIONS AUXILIARS ######
###############################


# Declaració de constants #

AA <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T",
        "V","W","Y")
PrPositions  <- paste0("P",1:99)
DRUGS  <-c("FPV","ATV","IDV","LPV", "NFV", "SQV" ,"TPV","DRV" )

# Càrrega de les posicions de resistència #
WTSeq <- read.table("Data/outputPI",
                    header=FALSE, col.names=paste0("P",1:99),sep="",colClasses = c("character"),
                    stringsAsFactors=TRUE)


# Declaració de funcions #

calc_freq <- function(x) return (lapply(x[1:ncol(x)],function(x)prop.table(table(x[!is.na(x)]))))


toSeq <- function(X) matrix(apply(X,1,function(x)paste(x,collapse="")),nrow=nrow(X),ncol=1)

toFactor <- function(X) {
  if(!methods::is(X,"data.frame"))  {
    X <- as.data.frame(X,stringsAsFactors = T)
  } else {
    X[] <- lapply(X,as.factor)
  }
  return(X)
}


pos_interac <- function(x) {
  comb2 <- t(combn(x,2))
  comb2 <- paste0(comb2[,1],comb2[,2])
  return(comb2)
}

rearrange_kint <- function(feat_space,rows,AA,PrPositions) {
  feat_space <- matrix(feat_space, nrow(feat_space), length(AA)*length(PrPositions))
  rownames(feat_space) <- rows
  colnames(feat_space) <- paste(rep(PrPositions,each=length(AA)),rep(AA,length(PrPositions)),sep="_")
  feat_space <- desparsify(feat_space,2)
  feat_space
}


#euclidean distance
calc_dist <- function(x,y) return(sqrt(sum((x - y)^2)))


vect_grep <- Vectorize(grepl)


searchSubs <- Vectorize(stringi::stri_count,vectorize.args = "fixed")

strSort <- function(x)  {  #Ordenar una string
  sapply(lapply(strsplit(x, NULL), sort), paste, collapse="")
}


permute_rep <- function(alphabet,l) { # Permutacions amb repetició
  topermute <- as.data.frame(matrix(rep(alphabet,l),nrow=length(alphabet),ncol=l))
  return(apply(expand.grid(topermute),1,function(x)paste(x,collapse = "")))
}


expand.grid.mod <- function(x, rep=FALSE) {
  g <- function(i) {
    z <- setdiff(x, x[seq_len(i-rep)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}


colourx <- function(x,WTSeq,DATAdum_delete) {
  wtalleles <-   combn( paste0("P",1:99,"_",WTSeq),2)
  wtalleles1 <-  sub(".*_", "", wtalleles)
  wtalleles2 <-  sub("_.*", "", wtalleles)
  wtalleles <- paste(paste(wtalleles2[1,],wtalleles2[2,],sep="_"), paste0(wtalleles1[1,],wtalleles1[2,]) ,sep="_")
  mutalleles <- colnames(DATAdum_delete)[!(colnames(DATAdum_delete) %in% paste0("P",1:99,"_",WTSeq))]
  mutalleles <- combn( mutalleles,2)
  mutalleles1 <-  sub(".*_", "", mutalleles)
  mutalleles2 <-  sub("_.*", "", mutalleles)
  mutalleles <- paste(paste(mutalleles2[1,],mutalleles2[2,],sep="_"), paste0(mutalleles1[1,],mutalleles1[2,]) ,sep="_")
  colors_wt <-  rep("grey90",length(x))
  colors_wt[rownames(x) %in% wtalleles] <- "lightskyblue"
  colors_wt[rownames(x)  %in% mutalleles] <- "wheat2"
  names(colors_wt) <- rownames(x)
  return(colors_wt)
}


plotldsign <- function(x,color) {
  plotImp(Rsquared[grep(x,names(Rsquared))]*D[grep(x,names(D))],main=x,
          nfeat=10,leftmargin = 7,color=color,absolute=FALSE,relative=FALSE)

}

