##########################
##### PREPROCESSING ######
##########################

# Declaració de funcions #

strSort <- function(x)    #Ordenar una string
  sapply(lapply(strsplit(x, NULL), sort), paste, collapse="")

strReverse <- function(x)   #Invertir una string
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

MixAASample <- function(x) {    # Mostreig de les barreges d'aminoàcids. Distribució uniforme.
  sample(unlist(strsplit(x,split="")),1,replace=TRUE)
}

is.trunc <- function(x) {
  return(grepl("\\*", x))
}

is.ins <- function(x) {
  return(grepl("\\#", x))
}

is.del <- function(x) {
  return(grepl("\\~", x))
}

is.amb <- function(x) {
  return(grepl("\\.|X|U|B|Z", x))
}



# Càrrega de dades #

DB <- c("PI", "NRTI","NNRTI", "INI")

for(NAME in DB) {
  print(NAME)

  # Per llegir directament la base de dades de la web
  # Raw <- read.table(paste0("http://hivdb.stanford.edu/download/GenoPhenoDatasets/",NAME,"_DataSet.Full.txt"),
  #                   header=TRUE, sep="\t", comment.char="@",  stringsAsFactors=FALSE)
  # Per llegir del fitxer descarregat:
  Raw <- read.table(paste0("Data/",NAME,"_DataSet.Full.txt"),
                     header=TRUE, sep="\t", comment.char="@",  stringsAsFactors=FALSE)


  if(NAME == "PI") {    # BASE DE DADES: PI (INHIBIDORS DE LA PROTEASA)
    Raw <- Raw[,-ncol(Raw)]
    PrPositions <- paste0("P",1:99)
    ### Consensus subtype B
    WTSeq <- read.table("Data/outputPI",
                        header=FALSE, col.names=names(Raw[,PrPositions]),sep="",colClasses = c("character"),
                        stringsAsFactors=TRUE)
  }

  if(NAME == "NRTI") {   # BASE DE DADES: NRTI (INHIBIDORS nucleòsids de la TRANSCRIPTASA REVERSA)
    Raw <- Raw[,-ncol(Raw)]
    PrPositions <- paste0("P",1:560)
    ### Consensus subtype B
    WTSeq <- read.table("Data/outputRT",
                        header=FALSE, col.names=names(Raw[,PrPositions]),sep="",colClasses = c("character"),
                        stringsAsFactors=TRUE)
  }

  if(NAME == "NNRTI")  { # BASE DE DADES: NNRTI (INHIBIDORS NO-nucleòsids de la TRANSCRIPTASA REVERSA)
    Raw <- Raw[,-ncol(Raw)]
    PrPositions <- paste0("P",1:560)
    ### Consensus subtype B
    WTSeq <- read.table("Data/outputRT",
                        header=FALSE, col.names=names(Raw[,PrPositions]),sep="",colClasses = c("character"),
                        stringsAsFactors=TRUE)
  }

  if(NAME == "INI")  { # BASE DE DADES: INI (INHIBIDORS INTEGRASA)
    Raw <- Raw[,-ncol(Raw)]
    PrPositions <- paste0("P",1:288)
    ### Consensus subtype B
    WTSeq <- read.table("Data/outputINI",
                        header=FALSE, col.names=names(Raw[,PrPositions]),sep="",colClasses = c("character"),
                        stringsAsFactors=TRUE)
  }


  # seqIDwt <- data.frame(0)
  # names(seqIDwt) <- "SeqID"
  # WTSeq <-cbind(seqIDwt,WTSeq)


  # Estandarditzar codis d'aminoàcids #

  ## Hi ha una problema previ amb els aminoàcids en barreja:
  ## Poden haver-hi diversos codis a la mateixa base de dades per a la mateixa barreja d'aa.
  ## Decidim ordenar alfabéticament les lletres de cada codi.
  ## Ex: IV i VI, per codificar una única mescla d'aminoàcids

  Raw.PrPositions <- Raw[,PrPositions]
  Raw.PrPositions[] <- apply(Raw.PrPositions,c(1,2),strSort)

  ## El mateix passa amb les lletres minúscules (ex: tenim certs casos de "l" pes compte de "L")

  Raw.PrPositions[] <- apply(Raw.PrPositions,c(1,2),toupper)

  ##  Per últim, pot ser recomanable canviar les "-" per la seqüència de la proteïna WT.

  for (j in 1:ncol(Raw.PrPositions)) {
    for (i in 1:nrow(Raw.PrPositions))
    {
      if (Raw.PrPositions[i,j] %in% c(".", "-")) Raw.PrPositions[i,j] <- WTSeq[j]
    }
  }
  Raw <- cbind(Raw[,setdiff(colnames(Raw),PrPositions)],Raw.PrPositions,stringsAsFactors=FALSE)

  write.csv(Raw, file = paste0("Data/",NAME,"_Seq.csv"),row.names = F)
  rm(Raw.PrPositions,Raw,i,j)
}

# Mostreig mescles d'aa #

# Funció principal de processament

processa.dades <- function(database,  mix.alleles = "random") {

  ## 1. Carrega la base de dades
  if(database=="PI"){
    PrPositions <- paste0("P",1:99)
    DRUGS  <-c("FPV","ATV","IDV","LPV", "NFV", "SQV" ,"TPV","DRV" )
  }
  if(database=="NRTI"){
    PrPositions <- paste0("P",1:560)
    DRUGS  <- c("X3TC","ABC","AZT","D4T", "DDI","TDF")
  }
  if(database=="NNRTI"){
    PrPositions <- paste0("P",1:560)
    DRUGS <- c("EFV" ,"NVP" ,"ETR", "RPV")
  }

  if(database=="INI"){
    PrPositions <- paste0("P",1:288)
    DRUGS  <- c("RAL","EVG","DTG","BIC","CAB")
  }

  Raw <- read.csv(paste0("Data/",database, "_Seq.csv"),colClasses=c("numeric","numeric","character","character","numeric",rep("character",3),
                                                   rep("numeric",length(DRUGS)),rep("character",length(PrPositions))))


  ## 2.Tractament de les posicions amb més d'un aa -> mostreig aleatori  o "res"

  if(mix.alleles=="random"){

    ## Mostreig aleatori: 1df

    Raw.notMix <- Raw[,PrPositions]
    Raw.notMix[] <- apply(Raw.notMix,c(1,2),MixAASample)
    Raw2 <- cbind(Raw[,setdiff(colnames(Raw),PrPositions)],Raw.notMix,stringsAsFactors=FALSE)
  }

  else {
    Raw2 <- Raw
  }

  ## 3. Eliminar ~, # i *

  delete.positions1 <-  which(apply(Raw[,PrPositions],c(1,2),is.trunc),arr.ind=TRUE) ## Moltes estan en mescles d'aa.
  delete.positions2 <-  which(apply(Raw[,PrPositions],c(1,2),is.ins),arr.ind=TRUE) ## No sabem quina és la inserció?
  delete.positions3 <-  which(apply(Raw[,PrPositions],c(1,2),is.del),arr.ind=TRUE)
  delete.positions4 <-  which(apply(Raw[,PrPositions],c(1,2),is.amb),arr.ind=TRUE)
  delete.positions <- c(delete.positions1[,1],delete.positions2[,1],delete.positions3[,1],delete.positions4[,1])
  if(length(delete.positions) !=0)  {
    delete.rows <- sort(unique(delete.positions))
    DATA <- Raw2[-delete.rows,]
  }
  else {
    DATA <- Raw2
  }

  # print(delete.positions)
  #
  # if(length(delete.positions) !=0)  {
  #   delete.rows <- unique(delete.positions[,1])
  #   Raw3 <- Raw2[-delete.rows,]
  # }
  # else {
  #   Raw3 <- Raw2
  # }


  # delete.positions <- which(Raw3=="#" | Raw3=="~", arr.ind = TRUE)
  # if(length(delete.positions) !=0)  {
  #     delete.rows <- unique(delete.positions[,1])
  #     DATA <- Raw3[-delete.rows,]
  # }
  # else {
  #     DATA <- Raw3
  # }

  ## 4. NAs -> esborrar
  ## Convertir · i X a NAs

  # DATA.PrPositions <- DATA[,PrPositions]
  # DATA.PrPositions[DATA.PrPositions=="." | DATA.PrPositions == "X" | DATA.PrPositions == "U" | DATA.PrPositions == "B"| DATA.PrPositions == "Z"] <- NA
  # print( which(DATA.PrPositions != "A" & DATA.PrPositions != "C" & DATA.PrPositions != "D" & DATA.PrPositions != "E" & DATA.PrPositions != "F" &
  #          DATA.PrPositions != "G" & DATA.PrPositions != "H" & DATA.PrPositions != "I" & DATA.PrPositions != "K" & DATA.PrPositions != "L" &
  #          DATA.PrPositions != "M" & DATA.PrPositions != "N" & DATA.PrPositions != "P" & DATA.PrPositions != "Q" & DATA.PrPositions != "R" &
  #          DATA.PrPositions != "S" & DATA.PrPositions != "T" & DATA.PrPositions != "V" & DATA.PrPositions != "W" & DATA.PrPositions != "Y",arr.ind = T) )## Revisió de les posicions
  #
  # DATA <- cbind(DATA[,setdiff(colnames(DATA),PrPositions)],DATA.PrPositions,stringsAsFactors=FALSE)
  # DATA <- DATA[complete.cases(DATA[,PrPositions]),]
  # print(dim(DATA))

  ## 4. Prendre log_10 de la variable target
  for(drug in DRUGS)  DATA[which(DATA[,drug] == 0),drug] <- min(DATA[,drug],na.rm = TRUE)/10
  DATA[,DRUGS] <- round(log10(DATA[,DRUGS]),digits = 6)
  for(drug in DRUGS) DATA[which(is.infinite(DATA[,drug] )),drug] <- NA

  ## 5.   Passar a factor
  # PrPositions.factor <-  lapply(DATA[,c("SeqID",PrPositions)], factor)
  # DATA <- merge(DATA[,setdiff(colnames(DATA),PrPositions)],PrPositions.factor,by="SeqID")
  DATA[,PrPositions] <-  lapply(DATA[, PrPositions], factor)

  ## 6. Guarda les dades processades
  return(DATA)

}


for(NAME in c("NRTI","NNRTI","INI")) {

  ## Versió sense mostrejar al·lels múltiples
  DATA <- processa.dades(database=NAME, mix.alleles = "no")
  save(list="DATA",file=paste0("Data/",NAME,  "_MultiA_clean"))

  ## Versió mostreig al·lels múltiples.
  DATA.list <- list()
  # bp.list <- list()
  for (i in 2:30) {
    DATA.l <- processa.dades(database=NAME, mix.alleles = "random")
    # PrPositions.factor <-  lapply(DATA.l[,c("SeqID",PrPositions)], factor)
    # DATA.l <- merge(DATA[,setdiff(colnames(DATA),PrPositions)],PrPositions.factor,by="SeqID")
    DATA.list[[i]] <- DATA.l
    # bp.list[[i]] <- badpos
    # print(dim(DATA.l))
  }

  save(list="DATA.list",file=paste0("Data/",NAME, "_Seq_clean"))
}
