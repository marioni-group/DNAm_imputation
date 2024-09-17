

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



targets = readRDS("GS20k_Targets_18869.rds")
agesex = read.csv("agemonths.csv")
pc = read.table("GS20K_ALL_MAF5_PCA.eigenvec")
dat <- readRDS("mvals.rds")

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
##########################################

lbc = readRDS('LBC_betas_3489_bloodonly.rds')

beta2m <- function(beta) {
  m <- log2(beta/(1-beta))
  return(m)
}

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

outlierID <- function(x, cut=4) {
  xx <- scale(x)
  retval <- ifelse (abs(xx) >= cut, TRUE, FALSE)
  retval
}

outlierTrim <- function(x, cut=4) {
  id <- outlierID(x, cut=cut)
  retval <- ifelse(id, NA, x)
}

cigs$pack_years_clean <- outlierTrim(log(cigs$pack_years + 1))
cigs <- cigs[,c(1:4,8:10)]

ahrr_lbc <- merge(ahrr_lbc, cigs, by.x="id", by.y="Basename")

cor.test(ahrr_lbc$ahrr_episcore, ahrr_lbc$ahrr_DNAm) 
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

