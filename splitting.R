###############################################
##### INCLUSION & EXCLUSION OF SEQUENCES ######
##### TRAINING AND TEST INDICES          ######
###############################################

# Càrrega de dades #

NAME <- "PI"
PrPositions <- paste0("P",1:99)

load(paste0("Data/",NAME,  "_MultiA_clean"))
rownames(DATA) <- 1:nrow(DATA)

##### Resum de les metadades
#### Nombre Seqüències per pacients
summary(factor(DATA$PtID))
# summary(factor(Raw$PtID[filt_seqs]))

#### Nombre seqüències per subtipus
(table(factor(DATA$Subtype)))
# (table(factor(Raw$Subtype[filt_seqs])))

#### Nombre seqüències per mètode de càlcul resistència
(table(factor(DATA$Method)))
# (table(factor(Raw$Method[filt_seqs])))

#### Nombre seqüències per Lab o Clinical
(table(factor(DATA$Type)))


### Filtrar: per subtipus i per metòde de calcular la resistència. Clinical perquè Lab són totes del mateix pacient:
work_idx <- (DATA$Subtype == "B") & (DATA$Method == "PhenoSense") & (DATA$Type == "Clinical")
nrow(DATA[work_idx,])
summary(factor(DATA[work_idx,"Subtype"]))
summary(factor(DATA[work_idx,"Method"]))

barplot(sort(table(factor(DATA$PtID[work_idx]))))


#### Seqüències amb ID repetits -> tenen RefID diferents
seqrep <-  names(which(table(factor(DATA$SeqID[work_idx])) >1))
View(subset(DATA,SeqID %in% seqrep))

View(subset(DATA,SeqID %in% seqrep)) ## la columna 5 és la de RefID
seqrep_idx <- rownames(subset(DATA,SeqID %in% seqrep)) ## la columna 5 és la de RefID
View(unique(DATA[seqrep_idx,-5]))
identical(sort(seqrep),sort(as.character(unique(DATA[seqrep_idx,-5])[,1])))

delete <- as.numeric(seqrep_idx[!(seqrep_idx %in% rownames(unique(DATA[seqrep_idx,-5])))])

work_idx[delete] <- FALSE

names(which(table(factor(DATA$SeqID[work_idx])) >1))
# character(0)
dim(DATA[work_idx,])

#### Pacients amb més d'1 seqüència:
patrep <- names(which(table(factor(DATA$PtID[work_idx])) >1))
View(subset(DATA[work_idx,],PtID %in% patrep))

patrep_dupl <- rep(0,length(patrep))
names(patrep_dupl) <- patrep
for(patient in patrep)  patrep_dupl[patient] <- anyDuplicated(subset(DATA[work_idx,],PtID == patient)[,PrPositions])
patrep_dupl <- patrep_dupl[patrep_dupl!=0]

#### Lo seu seria agafar les duplicades que tenen més fàrmacs!!
View(subset(DATA[work_idx,],PtID %in% names(patrep_dupl)))
delete <- c()
for(i in names(patrep_dupl)) {
  subset_data <- subset(DATA[work_idx,],PtID %in% i)
  revers <- duplicated(subset_data[,PrPositions],fromLast = TRUE)
  anvers <- duplicated(subset_data[,PrPositions],fromLast = FALSE)
  delete <- c(delete,rownames(subset_data[revers,,drop=F]))
  print(  subset_data[which(revers|anvers),])
  ## En general crec que em quedaria amb la darrera entrada, ja que normalment té més fàrmacs que la primera
}
delete <- as.numeric(delete)
work_idx[delete] <- FALSE

dim(DATA[work_idx,])
#  2087  115



# Test idx #

IDX <- list()
IDX$work <- which(work_idx)

## Sample: 20% sobre pacients
summary(factor(DATA[work_idx,"PtID"]))
sum(table(factor(DATA[work_idx,"PtID"]))>1)
# 243

patients <- unique(DATA[work_idx,"PtID"])
N <- length(patients)
Ntest <- round(0.20*N)
te_idx_pat <- sample(patients)[1:Ntest]

IDX$test <-  IDX$work[which(DATA[work_idx,"PtID"] %in% te_idx_pat)]

View(DATA[IDX$test,])



## Additional Test indexes ##

### HIV diferent del grup B
te_idx_subtype <- which(DATA$Subtype != "B" & DATA$Method =="PhenoSense")
length(te_idx_subtype)


### Antivirogram
te_idx_antivir <- which(DATA$Method =="Antivirogram" & DATA$Subtype == "B")
length(te_idx_antivir)

IDX$additional_test_subtype <- te_idx_subtype
IDX$additional_test_antivir <- te_idx_antivir
save(IDX,file=paste0("Partitions/",NAME,"_IDX"))
