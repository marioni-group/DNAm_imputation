

## load required libraries

library(tidyverse)
library(biglasso)

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
