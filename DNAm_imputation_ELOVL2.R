

## load required libraries

library(bigmemory)
library(bigmemoryExt)
library(optparse)
library(biglasso)
library(tidyverse)
library(foreign)
library(coxme)
library(kinship2)
library(preputils)


## define functions

cbindBM_list <- function(x, binding="right", 
                         z=NULL, type=NULL, separated=NULL,
                         backingfile=NULL, backingpath=NULL,
                         descriptorfile=NULL, binarydescriptor=FALSE,
                         shared=TRUE, erase = TRUE)
{
  
  if (is.null(type)) type <- typeof(x[[1]])
  if (is.big.matrix(x[[1]])) {
    if (is.null(separated)) separated <- is.separated(x[[1]])
  } else {
    separated <- FALSE
  }
  
  cols_list <- list()
  total_cols <- 0
  for (i in 1:length(x)) {
    cols <- cleanupcols(NULL, ncol(x[[i]]), colnames(x[[i]]))
    cols_list <- append(cols_list, list(cols))
    total_cols <- total_cols + ncol(x[[i]])
  }    
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(x[[1]]), ncol=total_cols, type=type, init=NULL,
                    dimnames=dimnames(x[[1]]), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared=shared)
  }
  
  counter <- 0
  for (i in 1:length(cols_list)) {
    print(i)
    if (i == 1) {
      z[, 1:length(cols_list[[i]])] <- x[[i]][,cols_list[[i]]]
    } else {
      z[, (counter + 1):(counter + length(cols_list[[i]]))] <- x[[i]][,cols_list[[i]]]
    }
    counter <- counter + length(cols_list[[i]])
    print(counter)
    
    if (erase == TRUE) {
      cat("\nErasing chunk and liberating memory...\n\n")
      x[[i]] <- "Replacement"
      gc()
    }
  }
  return(z)
}

cleanupcols <- function(cols=NULL, nc=NULL, colnames=NULL) {
  if (is.null(cols)) cols <- 1:nc
  else {
    if (!is.numeric(cols) & !is.character(cols) & !is.logical(cols))
      stop("column indices must be numeric, logical, or character vectors.")
    if (is.character(cols))
      if (is.null(colnames)) stop("column names do not exist.")
    else cols <- mmap(cols, colnames)
    if (is.logical(cols)) {
      if (length(cols) != nc)
        stop(paste("column vector length must match the number of",
                   "columns of the matrix."))
      cols <- which(cols)
    }
    tempj <- .Call("CCleanIndices", as.double(cols), as.double(nc), PACKAGE="bigmemory")
    if (is.null(tempj[[1]])) stop("Illegal column index usage in extraction.\n")
    if (tempj[[1]]) cols <- tempj[[2]]
  }
  return(cols)
}



targets = readRDS("GS20k_Targets_18869.rds")
agesex = read.csv("agemonths.csv")
pc = read.table("GS20K_ALL_MAF5_PCA.eigenvec")
dat <- readRDS("mvals.rds")

elovl2_cpg <- dat[which(rownames(dat)=="cg16867657"),]
agesex = agesex[which(agesex$id %in% targets$Sample_Name),]
elovl2_cpg <- as.data.frame(elovl2_cpg)
elovl2_cpg$id <- rownames(elovl2_cpg)
elovl2_cpg = elovl2_cpg[which(elovl2_cpg$id %in% targets$X), c("id", "elovl2_cpg")]
ids <- targets[,1:2]
elovl2_cpg <- merge(elovl2_cpg, ids, by.x="id", by.y="X", all=TRUE)
names(elovl2_cpg) <- c("meth_id", "elovl2_cpg", "id")
pheno = agesex
pheno$elovl2_cpg = elovl2_cpg[match(agesex$id, elovl2_cpg$id), "elovl2_cpg"]
pheno = merge(pheno, pc, by.x="id", by.y="V2", all.x=TRUE)
batch <- targets[,c("Sample_Name", "Batch")]
pheno <- merge(pheno, batch, by.x="id", by.y="Sample_Name", all=TRUE)
agesex <- read.table("../../agesex.txt", header=TRUE)
pheno$elovl2_cpg[!(pheno$id %in% agesex$id)] <- NA
pheno$age_months[!(pheno$id %in% agesex$id)] <- NA
pheno$elovl2_cpg_res = resid(lm(elovl2_cpg ~ sex + as.factor(Batch) + V3 + V4 + V5 + V6 +
                     V7 + V8 + V9 + V10 + V11 + V12, data=pheno, na.action=na.exclude))

####################################
### PREP DNAm 
####################################

anno <- readRDS("annotations.rds")

dat <- dat[-which(grepl("cg16867657", rownames(dat))),]

# NA's - mean-impute data
for(i in 1:ncol(dat)){
    if(length(which(is.infinite(dat[,i])))>=1){
    dat[which(is.infinte(dat[,i])),i] = NA
  }

  if(length(which(is.na(dat[,i])))>=1){
    dat[which(is.na(dat[,i])),i] = mean(dat[,i], na.omit=T)
  }
}

# filter probes to those on 450K array
probes_450 <- rownames(anno[which(anno$Methyl450_Loci==TRUE),])

dat <- dat[which(rownames(dat) %in% probes_450),]

# find top 200k most variable probes
var_probes = apply(dat, 1, sd)
probes_200k = names(var_probes)[rev(order(var_probes))]

cpg = rownames(dat)

pheno$Sample_Sentrix_ID = targets[match(pheno$id, targets$Sample_Name), "Sample_Sentrix_ID"]
dat = dat[which(rownames(dat) %in% cpg),pheno$Sample_Sentrix_ID]

# filter probes to 200k variable probes
probes_200k_450 = probes_200k[which(probes_200k %in% rownames(dat))][1:200000]
dat = dat[probes_200k_450, ]


elovl2_cpg <- pheno
elovl2_cpg = elovl2_cpg[-which(is.na(elovl2_cpg$elovl2_cpg_res)), ]
meth <- dat[,elovl2_cpg$Sample_Sentrix_ID]
div <- 15 # Number of chunks to divide OG methylation dataframe
por <- ceiling(length(colnames(meth))/div)
chunk_list <- list()
for (i in 1:div) {
  cat(paste0("\nWorking on chunk: ", i, " of ", div))
  if (i == 1) {
    chunk <- as.big.matrix(meth[,1:(por-1)])
  } else if (i == div) {
    chunk <- as.big.matrix(meth[,(por*(i-1)):length(colnames(meth))])
  } else {
    chunk <- as.big.matrix(meth[,(por*(i-1)):((por*i)-1)])
  }
  cat("\nMade chunk. Appending to chunk list...\n")
  chunk_list <- append(chunk_list, list(chunk))
  gc()
}
# Saving names prior to chunk fusing
names <- colnames(meth)
rm(meth)
cat("\nRAM clean up...\n\n")
gc()
cat("\nFusing chunks!\n\n")
dat2 <- cbindBM_list(x = chunk_list)
rm(chunk, chunk_list)
# Set CpG names
options(bigmemory.allow.dimnames=TRUE)
colnames(dat2)<- names
dat2 = t(dat2)
colnames(dat2)= rownames(dat)
rownames(dat2) <- names

####################################
### RUN TRAINING
####################################

location <- "methylation_lasso/"
set.seed(1783) # set seed to ensure fold variation minimised
# Set list of traits to run through
list <- "elovl2_cpg_residual"

elovl2_cpg$Sample_Sentrix_ID = targets[match(elovl2_cpg$id, targets$Sample_Name), "Sample_Sentrix_ID"]
elovl2_cpg = elovl2_cpg[match(rownames(dat2), elovl2_cpg$Sample_Sentrix_ID),]
identical(rownames(dat2), elovl2_cpg$Sample_Sentrix_ID) # TRUE
ytrain <- elovl2_cpg
i <- "elovl2_cpg_res"
q <- ytrain[,"Sample_Sentrix_ID"] # Get just basenames for people in the y variable
p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable
name_p <- colnames(p) # Get the name of the trait for this iteration
y <- cbind(q,p) # Bind Basename and protein data together into one set
names(y)[2] <- "pheno" # Assign a generic name to the protein variable
y <- as.numeric(y$pheno) # Create a numeric list for this variable to feed in as y to the model
# Run training
lasso.cv <- cv.biglasso(dat2, y, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
fit <- biglasso(dat2, y, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lasso.cv$lambda.min) # model fit
elovl2_coefs <- coef(fit) # Extract coeficients
elovl2_coefs <- as.data.frame(elovl2_coefs[which(elovl2_coefs!=0),]) # Remove coeficients that are 0
elovl2_coefs$Predictor <- name_p # Assign protein identifier
names(elovl2_coefs)[1] <- "Coefficient" # Tidy naming
elovl2_coefs$CpG <- rownames(elovl2_coefs) # Create episcores column
write.csv(elovl2_coefs, file = paste0(location, name_p, "_GS20k_train_weights_cg16867657_200k.csv"), row.names = F)

elovl2_coefs$chr = anno[match(elovl2_coefs$CpG, anno$Name), 'chr']
elovl2_coefs$pos = anno[match(elovl2_coefs$CpG, anno$Name), 'pos']
table(elovl2_coefs$chr)
head(elovl2_coefs[order(abs(elovl2_coefs$Coefficient), decreasing=T),])


##########################################
### TEST COEFFICIENTS IN EXTERNAL DATASET
### Hannum
##########################################

hannum <- readRDS("methylation.RDS")

beta2m <- function(beta) {
  m <- log2(beta/(1-beta))
  return(m)
}

hannum[hannum==1] = 0.999999999
hannum[hannum==0] = 0.000000001
hannum <- t(hannum)

elovl2_hannum = beta2m(hannum['cg16867657',])

probes = elovl2_coefs$CpG
probes = probes[-1]

# extract coefficient probes from dataset
hannum = hannum[rownames(hannum) %in% probes,]
hannum = beta2m(hannum)

elovl2_coefs_hannum = elovl2_coefs[match(rownames(hannum), elovl2_coefs$CpG),]
all.equal(rownames(elovl2_coefs_hannum), rownames(hannum))

hannum_elovl2 = data.frame(id = colnames(hannum))

# calculate episcores
score = list()
for(i in 1:nrow(hannum_elovl2)){
  score[[i]] = sum(hannum[,i] * elovl2_coefs_hannum$Coefficient) + elovl2_coefs[grep('Intercept', elovl2_coefs$CpG),'Coefficient']
}

hannum_elovl2$elovl2_episcore = unlist(score)

elovl2_hannum <- as.data.frame(elovl2_hannum)
names(elovl2_hannum) <- "elovl2_DNAm"
elovl2_hannum$id <- rownames(elovl2_hannum)
elovl2_hannum <- merge(elovl2_hannum, hannum_elovl2, by="id")

## test correlation of episcores with age
hannum_age <- readRDS("samples.RDS")

hannum_age <- hannum_age[,c("geo_accession", "age (y):ch1")]
names(hannum_age)[2] <- "age"

elovl2_hannum <- merge(elovl2_hannum, hannum_age, by.x="id", by.y="geo_accession", all.x=TRUE)
elovl2_hannum$age <- as.numeric(elovl2_hannum$age)

cor.test(elovl2_hannum$elovl2_episcore, elovl2_hannum$elovl2_DNAm) 
cor.test(elovl2_hannum$elovl2_episcore, elovl2_hannum$age)
cor.test(elovl2_hannum$elovl2_DNAm, elovl2_hannum$age)

write.csv(elovl2_hannum, "elovl2_cg16867657_Hannum_measured_and_predicted.csv", row.names=FALSE, quote=FALSE)


##########################################
### TEST COEFFICIENTS IN EXTERNAL DATASET
### GSE246337
##########################################

boac <- as.data.frame(fread("zcat GSE246337_betas.csv.gz", header=TRUE))

# remove suffix from probe names
boac$V1 <- gsub("_.*", "", boac$V1)

elovl2_boac_dup = boac[which(boac$V1=="cg16867657"),]

# remove duplicate probes and set rownames as probe names
boac <- boac[!duplicated(boac$V1),]
rownames(boac) <- boac$V1
boac <- boac[,2:501]
boac <- as.matrix(boac)

elovl2_boac = beta2m(boac['cg16867657',])

probes = elovl2_coefs$CpG
probes = probes[-1]

# extract coefficient probes from dataset
boac = boac[rownames(boac) %in% probes,]
boac = beta2m(boac)

elovl2_coefs_boac = elovl2_coefs[match(rownames(boac), elovl2_coefs$CpG),]
all.equal(rownames(elovl2_coefs_boac), rownames(boac))

boac_elovl2 = data.frame(id = colnames(boac))

# calculate episcores
score = list()
for(i in 1:nrow(boac_elovl2)){
  score[[i]] = sum(boac[,i] * elovl2_coefs_boac$Coefficient) + elovl2_coefs[grep('Intercept', elovl2_coefs$CpG),'Coefficient']
}

boac_elovl2$elovl2_episcore = unlist(score)

elovl2_boac <- as.data.frame(elovl2_boac)
names(elovl2_boac) <- "elovl2_DNAm"
elovl2_boac$id <- rownames(elovl2_boac)
elovl2_boac <- merge(elovl2_boac, boac_elovl2, by="id")

## test correlation of episcores with age
boac_age <- read.table("ages.txt", header=TRUE)

names(boac_age) <- c("id", "age")

elovl2_boac <- merge(elovl2_boac, boac_age, by="id", all.x=TRUE)
elovl2_boac$age <- as.numeric(elovl2_boac$age)

cor.test(elovl2_boac$elovl2_episcore, elovl2_boac$elovl2_DNAm) 
cor.test(elovl2_boac$elovl2_episcore, elovl2_boac$age)
cor.test(elovl2_boac$elovl2_DNAm, elovl2_boac$age)

write.csv(elovl2_boac, "elovl2_cg16867657_GSE246337_measured_and_predicted.csv", row.names=FALSE, quote=FALSE)
