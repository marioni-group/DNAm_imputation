

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
library(foreign)
library(stringr)
library(imputeTS)
library(data.table)
library(limma)
library(foreach)
library(doParallel)
library(caret)
library(pROC)


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

beta2m <- function(beta) {
  m <- log2(beta/(1-beta))
  return(m)
}

outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}



targets = readRDS("GS20k_Targets_18869.rds")
pc = read.table("GS20K_ALL_MAF5_PCA.eigenvec")

####################################
### AHRR
####################################

dat <- readRDS("mvals.rds")
agesex = read.csv("agemonths.csv")

agesex = agesex[which(agesex$id %in% targets$Sample_Name),]
ahrr_cpg <- dat[which(rownames(dat)=="cg05575921"),]
ahrr_cpg <- as.data.frame(ahrr_cpg)
ahrr_cpg$id <- rownames(ahrr_cpg)
ahrr_cpg = ahrr_cpg[which(ahrr_cpg$id %in% targets$X), c("id", "ahrr_cpg")]
ids <- targets[,1:2]
ahrr_cpg <- merge(ahrr_cpg, ids, by.x="id", by.y="X", all=TRUE)
names(ahrr_cpg) <- c("meth_id", "cpg", "id")
pheno = agesex
pheno$ahrr_cpg = ahrr_cpg[match(agesex$id, ahrr_cpg$id), "cpg"]
pheno = merge(pheno, pc, by.x="id", by.y="V2", all.x=TRUE)
batch <- targets[,c("Sample_Name", "Batch")]
pheno <- merge(pheno, batch, by.x="id", by.y="Sample_Name", all=TRUE)
agesex <- read.table("agesex.txt", header=TRUE)
pheno$ahrr_cpg[!(pheno$id %in% agesex$id)] <- NA
pheno$age_months[!(pheno$id %in% agesex$id)] <- NA
pheno$ahrr_cpg_res = resid(lm(ahrr_cpg ~ sex + as.factor(Batch) + V3 + V4 + V5 + V6 +
                     V7 + V8 + V9 + V10 + V11 + V12, data=pheno, na.action=na.exclude))

####################################
### PREP DNAm 
####################################

anno <- readRDS("annotations.rds")

dat <- dat[-which(grepl("cg05575921", rownames(dat))),]

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


ahrr_cpg <- pheno
ahrr_cpg = ahrr_cpg[-which(is.na(ahrr_cpg$ahrr_cpg_res)), ]
meth <- dat[,ahrr_cpg$Sample_Sentrix_ID]
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
list <- "ahrr_cpg_residual"

ahrr_cpg$Sample_Sentrix_ID = targets[match(ahrr_cpg$id, targets$Sample_Name), "Sample_Sentrix_ID"]
ahrr_cpg = ahrr_cpg[match(rownames(dat2), ahrr_cpg$Sample_Sentrix_ID),]
identical(rownames(dat2), ahrr_cpg$Sample_Sentrix_ID) # TRUE
ytrain <- ahrr_cpg
i <- "ahrr_cpg_res"
q <- ytrain[,"Sample_Sentrix_ID"] # Get just basenames for people in the y variable
p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable
name_p <- colnames(p) # Get the name of the trait for this iteration
y <- cbind(q,p) # Bind Basename and protein data together into one set
names(y)[2] <- "pheno" # Assign a generic name to the protein variable
y <- as.numeric(y$pheno) # Create a numeric list for this variable to feed in as y to the model
# Run training
lasso.cv <- cv.biglasso(dat2, y, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
fit <- biglasso(dat2, y, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lasso.cv$lambda.min) # model fit
ahrr_coefs <- coef(fit) # Extract coeficients
ahrr_coefs <- as.data.frame(ahrr_coefs[which(ahrr_coefs!=0),]) # Remove coeficients that are 0
ahrr_coefs$Predictor <- name_p # Assign protein identifier
names(ahrr_coefs)[1] <- "Coefficient" # Tidy naming
ahrr_coefs$CpG <- rownames(ahrr_coefs) # Create episcores column
write.csv(ahrr_coefs, file = paste0(location, name_p, "_GS20k_train_weights_cg05575921_200k.csv"), row.names = F)

ahrr_coefs$chr = anno[match(ahrr_coefs$CpG, anno$Name), 'chr']
ahrr_coefs$pos = anno[match(ahrr_coefs$CpG, anno$Name), 'pos']
table(ahrr_coefs$chr)
head(ahrr_coefs[order(abs(ahrr_coefs$Coefficient), decreasing=T),] )


##########################################
### TEST COEFFICIENTS IN EXTERNAL DATASET
### LBC
##########################################

lbc = readRDS('LBC_betas_3489_bloodonly.rds')

ahrr_lbc = beta2m(lbc['cg05575921',])

probes = ahrr_coefs$CpG
probes = probes[-1]

# extract coefficient probes from dataset
lbc = lbc[rownames(lbc) %in% probes,]
lbc = beta2m(lbc)

ahrr_coefs_lbc = ahrr_coefs[match(rownames(lbc), ahrr_coefs$CpG),]
all.equal(rownames(ahrr_coefs_lbc), rownames(lbc))

lbc_ahrr = data.frame(id = colnames(lbc))

# calculate episcores
score = list()
for(i in 1:nrow(lbc_ahrr)){
  score[[i]] = sum(lbc[,i] * ahrr_coefs_lbc$Coefficient) + ahrr_coefs[grep('Intercept', ahrr_coefs$CpG),'Coefficient']
}

lbc_ahrr$ahrr_episcore = unlist(score)

ahrr_lbc <- as.data.frame(ahrr_lbc)
ahrr_lbc$id <- rownames(ahrr_lbc)
names(ahrr_lbc)[1] <- "ahrr_DNAm"
ahrr_lbc <- merge(ahrr_lbc, lbc_ahrr, by="id")


## test correlation of episcores with smoking phenotypes 
cigs <- read.spss("LBC1936_EWAS_BldBasedChonicLowGradeInflammation_RH_27FEB2023.sav",
    to.data.frame = TRUE,
    use.value.labels = FALSE)

lbc_targets <- readRDS("targets_3489_bloodonly.rds")

lbc_targets <- lbc_targets[which(lbc_targets$WAVE==1 & lbc_targets$cohort=="LBC36"),]

cigs <- merge(cigs, lbc_targets, by.x="lbc36no", by.y="ID_raw")
cigs <- cigs[,c("Basename", "lbc36no", "age", "sex.y", "smokagestart_w1", "smokagestop_w1", "smoknumcigs_w1")]
# create smoking status and pack years variables
cigs$smoke <- NA
cigs$smoke[is.na(cigs$smokagestart_w1) & is.na(cigs$smokagestop_w1)] <- 0 #never
cigs$smoke[!is.na(cigs$smokagestart_w1) & !is.na(cigs$smokagestop_w1)] <- 1 #former
cigs$smoke[!is.na(cigs$smokagestart_w1) & is.na(cigs$smokagestop_w1)] <- 2 #current
cigs$pack_years <- NA

for(i in 1:nrow(cigs)){
  if(cigs$smoke[i]==2){
    cigs$pack_years[i] <- (cigs$age[i] - cigs$smokagestart_w1[i]) * (cigs$smoknumcigs_w1[i] / 20)
  } else if(cigs$smoke[i]==1){
    cigs$pack_years[i] <- (cigs$smokagestop_w1[i] - cigs$smokagestart_w1[i]) * (cigs$smoknumcigs_w1[i] / 20)
  } else if(cigs$smoke[i]==0){
    cigs$pack_years[i] <- 0
  }
}

cigs$pack_years_clean <- outlierTrim(log(cigs$pack_years + 1))
cigs <- cigs[,c(1:4,8:10)]

ahrr_lbc <- merge(ahrr_lbc, cigs, by.x="id", by.y="Basename")

cor.test(ahrr_lbc$ahrr_episcore, ahrr_lbc$ahrr_DNAm) 
cor.test(ahrr_lbc$ahrr_episcore, ahrr_lbc$pack_years)
cor.test(ahrr_lbc$ahrr_DNAm, ahrr_lbc$pack_years)
cor.test(ahrr_lbc$ahrr_episcore, ahrr_lbc$pack_years_clean)
cor.test(ahrr_lbc$ahrr_DNAm, ahrr_lbc$pack_years_clean)

write.csv(ahrr_lbc, "ahrr_cg05575921_LBC_measured_and_predicted.csv", row.names=FALSE, quote=FALSE)

# calculate area under curve scores for smoking status
current_never <- subset(ahrr_lbc, smoke %in% c(0, 2))
former_never <- subset(ahrr_lbc, smoke %in% c(0, 1))
current_former <- subset(ahrr_lbc, smoke %in% c(1, 2))

current_never$current <- ifelse(current_never$smoke == 2, 1, 0)
former_never$former <- former_never$smoke
current_former$current <- ifelse(current_former$smoke == 2, 1, 0)

# episcore
roc36_current_never <- roc(response = current_never$current, predictor = current_never$ahrr_episcore)
auc_36_current_never <- auc(roc36_current_never)

roc36_former_never <- roc(response = former_never$former, predictor = former_never$ahrr_episcore)
auc_36_former_never <- auc(roc36_former_never)

roc36_current_former <- roc(response = current_former$current, predictor = current_former$ahrr_episcore)
auc_36_current_former <- auc(roc36_current_former)

# measured probe
roc36_current_never <- roc(response = current_never$current, predictor = current_never$ahrr_DNAm)
auc_36_current_never <- auc(roc36_current_never)

roc36_former_never <- roc(response = former_never$former, predictor = former_never$ahrr_DNAm)
auc_36_former_never <- auc(roc36_former_never)

roc36_current_former <- roc(response = current_former$current, predictor = current_former$ahrr_DNAm)
auc_36_current_former <- auc(roc36_current_former)




####################################
### ELOVL2
####################################

dat <- readRDS("mvals.rds")
agesex = read.csv("agemonths.csv")

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
boac_age <- read.table("GSE246337/ages.txt", header=TRUE)

names(boac_age) <- c("id", "age")

elovl2_boac <- merge(elovl2_boac, boac_age, by="id", all.x=TRUE)
elovl2_boac$age <- as.numeric(elovl2_boac$age)

cor.test(elovl2_boac$elovl2_episcore, elovl2_boac$elovl2_DNAm) 
cor.test(elovl2_boac$elovl2_episcore, elovl2_boac$age)
cor.test(elovl2_boac$elovl2_DNAm, elovl2_boac$age)

write.csv(elovl2_boac, "GSE246337/elovl2_cg16867657_GSE246337_measured_and_predicted.csv", row.names=FALSE, quote=FALSE)



####################################
### GDF15
####################################

#### PREPARE THE DATA ####

opt_out <- read.csv('participant_opt_out.csv', header = FALSE)

# read in proteins and remove opt-outs
prots_long <- read.csv('olink_data_big.txt', sep = '\t')  %>%
  filter(ins_index == 0 & !eid %in% opt_out$V1) # select just baseline assessment

# protein_names (data coding 143; https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=143&nl=1)
prot_codings <- read.csv('coding143.tsv', sep = '\t')
# just keep name before semicolon
prot_codings$meaning <- sub(';.*', '', prot_codings$meaning)
prot_codings <- prot_codings %>%
  rename(protein_id = coding)

# add protein names instead of codings
prots_long <- merge(prots_long, prot_codings,
                    by = 'protein_id',
                    all.y = TRUE)
# to wide format
prots <- prots_long %>%
  select(eid, meaning, result) %>%
  pivot_wider(names_from = meaning,
              values_from = result) %>%
    rename(id = eid)

# read in the proteomics meta-data
prots_meta <- read.csv('proteomics_meta.csv') %>%
  select(-ends_with(c('.1.0', '.2.0', '.3.0'))) %>%
  rename(id = eid,
         protein_n = X30900.0.0,
         plate = X30901.0.0,
         well = X30902.0.0,
         ukb_ppp = X30903.0.0) %>%
  filter(!is.na(protein_n))
prots_meta$ukb_ppp[is.na(prots_meta$ukb_ppp)] <- 0

# merge with protein file
prots <- merge(prots_meta, prots, by = 'id', all.y = TRUE)

# check numbers of missing values for each column
missing_vals <- sapply(prots, function(x) sum(is.na(x)))
# find proteins missing in more than 20% of participants and remove them
bad_prots <- missing_vals[missing_vals > nrow(prots)/5]
prots_clean <- prots %>%
  select(-names(bad_prots)) %>%
  #filter(ukb_ppp == 0) %>% # potentially remove participants enriched for diseases
  select(-ukb_ppp) %>%
# find participants missing 20% or more proteins and remove them
  filter(protein_n/nrow(prot_codings) >= 0.8) %>%
  # remove those missing GDF15 data
  filter(!is.na(GDF15)) %>%
  select(-protein_n)

prots_clean$id <- as.factor(prots_clean$id)

# columns with dashes in name cause problems;
# we will save them into a vector and replace the dashes with underscores
columns_with_dash <- grep('-', names(prots_clean), value = TRUE)
columns_with_dash_list <- as.list(columns_with_dash)
new_column_names <- gsub('-', '_', columns_with_dash)
names(prots_clean)[names(prots_clean) %in% columns_with_dash] <- new_column_names

# split data into training and test sets
set.seed(24)
train_indices <- sample(1:nrow(prots_clean), round(nrow(prots_clean)/2))
train_data <- prots_clean[train_indices, ]
test_data <- prots_clean[-train_indices, ]

# impute missing values; separately for training and test data
train_imputed <- data.frame(train_data)
test_imputed <- data.frame(test_data)
relevant_columns <- colnames(select(prots_clean, -id, -plate, -well))

#for (col in relevant_columns) {
#  train_imputed[[col]][is.na(train_imputed[[col]])] <-
#   sample(train_imputed[[col]][!is.na(train_imputed[[col]])], 1)
#
#  test_imputed[[col]][is.na(test_imputed[[col]])] <-
#    sample(test_imputed[[col]][!is.na(test_imputed[[col]])], 1)
#}

# mean-imputation alternative
for (col in relevant_columns) {
  train_imputed[[col]][is.na(train_imputed[[col]])] <-
    mean(train_imputed[[col]], na.rm = TRUE)

  test_imputed[[col]][is.na(test_imputed[[col]])] <-
    mean(test_imputed[[col]], na.rm = TRUE)
}

# scale the predictors and the outcome
numeric_variables <- colnames(prots_clean[, c(sapply(prots_clean[, ], function(x) is.numeric(x)))])
train_imputed[numeric_variables] <- 
  sapply(train_imputed[, numeric_variables], function(x) as.vector(scale(x)))
test_imputed[numeric_variables] <- 
  sapply(test_imputed[, numeric_variables], function(x) as.vector(scale(x)))



#### TRAIN THE MODEL ####

library(biglasso)
y <- train_imputed[,which(names(train_imputed)=="GDF15")]
x1 <- as.matrix(train_imputed[,-c(1:3,which(names(train_imputed)=="GDF15"))])
X.bm <- as.big.matrix(x1)

lasso.cv <- cv.biglasso(X.bm, y, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
fit <- biglasso(X.bm, y, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lasso.cv$lambda.min) # model fit
gdf_coefs <- coef(fit) # Extract coeficients
gdf_coefs <- as.data.frame(gdf_coefs[which(gdf_coefs!=0),]) # Remove coeficients that are 0
names(gdf_coefs)[1] <- "Coefficient" # Tidy naming
gdf_coefs$prot <- rownames(gdf_coefs) # Create episcores column

location <- "Olink/"
write.csv(gdf_coefs[,c(2,1)], file = paste0(location, "_UKB_train_weights_gdf.csv"), row.names = F) 
  


#### PREPARE FOR MODELLING IN TEST DATA ####

cols <- which(names(test_imputed) %in% gdf_coefs$prot)
test <- test_imputed[,cols]

wts <- gdf_coefs$Coefficient[-1]

pred <- as.matrix(test) %*% wts

test_imputed$gdf_pred <- pred +  gdf_coefs$Coefficient[1]

# load the causes of censoring and merge them
dementia <- read.csv('ADOs_dementia.csv') %>%
  mutate(dementia_date = as.Date(X42018.0.0, format = '%Y-%m-%d')) %>%
  select(eid, dementia_date)
death <- read.csv('death_date.csv') %>%
  mutate(death_date = as.Date(X40000.0.0, format = '%Y-%m-%d')) %>%
  select(eid, death_date)
loss <- read.csv('loss_to_follow_up.csv') %>%
  mutate(loss_date = as.Date(X191.0.0, format = '%Y-%m-%d')) %>%
  select(eid, loss_date)

# determine censoring date
cens <- merge(dementia, death, by = 'eid')
cens <- merge(cens, loss, by = 'eid') %>%
  rename(id = eid) %>%
  transform(cens_date = pmin(dementia_date, death_date, loss_date, na.rm = TRUE))
cens$cens_date[is.na(cens$cens_date)] <- as.Date('2022-12-31', format = '%Y-%m-%d')

# create dementia outcome variable
cens$dementia <- 0
cens$dementia[!is.na(cens$dementia_date)] <- 1

# calculate follow-up
dates <- read.csv('dems.csv') %>%
  mutate(ass_date = as.Date(X53.0.0, format = '%Y-%m-%d')) %>%
  select(eid, ass_date) %>%
  rename(id = eid) %>%
  merge(., cens, by = 'id')
dates$follow_up <- as.numeric(difftime(dates$cens_date, dates$ass_date))/365.25

# merge with GDF15 data frame  
df <- merge(test_imputed[,c("id","gdf_pred","GDF15")], dates, by = 'id', all.x = TRUE)

# add age and sex
agesex <- read.csv('age_sex_formatted.csv', sep = '|') %>%
  select(id, sex, birth_date)
df <- merge(df, agesex, by = 'id', all.x = TRUE)
df$age <- as.numeric(difftime(df$ass_date, df$birth_date))/365.25

# remove negative follow-up (dementia before assessment),
# and people younger than 60 at censoring
df <- df %>%
  filter(follow_up > 0 & 
           birth_date < as.Date('1962-01-01', format = '%Y-%m-%d'))



#### RUN COX MODEL ####

# run Cox model
fit_hr_ref <- survival::coxph(data = df,
                          survival::Surv(follow_up, dementia) ~ scale(GDF15) + age + sex)
fit_hr_score <- survival::coxph(data = df,
                          survival::Surv(follow_up, dementia) ~ scale(gdf_pred) + age + sex)

broom::tidy(fit_hr_ref, exponentiate = T, conf.int = T)[1,]
broom::tidy(fit_hr_score, exponentiate = T, conf.int = T)[1,]
cor.test(df$GDF15, df$gdf_pred)

plot(df$GDF15, df$gdf_pred)
