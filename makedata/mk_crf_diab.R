################################################################################
# PROD PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/makedata/production/mk_crf_diab.r
# WORK PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/makedata/working/mk_crf_diab.r
# 
# RISK ASSESSMENT: Low
# DOCUMENT ID:     PROG.10084.v3
# 
# PURPOSE:         Predict individual treatment effects (ITEs) for Diabetes at yr 5
#                  using causal random forests (continuous risk factors)
# SOURCE PRGM:     none
# INPUT:           rdata/adgrf.rds
# OUTPUT:          rdata/crf_diab.rds (fit crf object)
#                  rdata/ite_diab.rds (predicted ites)
# PACKAGES USED:   tidyverse, grf, lmtest, sandwich, survival, caret, this.path, skimr
# FUNCTIONS USED:  write_contents
# EXEMPTIONS:      none
# 
# AUTHOR:          Emma Davies Smith
# CREATION DATE:   09FEB2026
# 
# NOTES:           
# MODIFICATIONS:    
################################################################################


################################################################################
# SETUP
################################################################################

# Project setupfile
source("/home/reprieve/POOLED/final/MACE_hte/programs/setupfile.R")

# Libraries
library(tidyverse)              # data manipulation
library(grf,       lib=pkgdir)  # generalized random forests
library(lmtest,    lib=pkgdir)  # dependency of grf that doesnt load properly
library(sandwich,  lib=pkgdir)  # dependency of grf that doesnt load properly 
library(survival)               # survival methods (default package)
library(caret,     lib=pkgdir)  # create dummy variables
library(this.path ,lib=pkgdir)  # file paths

# RWG: write_contents (requires skimr library)
library(skimr, lib=pkgdir)
source("/home/reprieve/POOLED/final/MACE_hte/functions/write_contents.R")


################################################################################
# LOAD DATA 
################################################################################

# Load dataset: adgrf
adgrf <- readRDS(paste0(rdata, 'adgrf.RDS'))

# Subset to complete cases 
adgrfcc <- adgrf %>% filter(ccfl) 

# Specify covariates of interest
# !!! NOTE: MAKE SURE THIS ALIGNS WITH covar SPECIFIED IN mk_adgrf.R !!!
# NOTE: Treatment arm is not included as part of covar, see below.
covar <- c('gbdrace', 'natalsex', 'age', 
           'sysbp', 'schol', 'shdl', 'bmi', 'egfr_eq21', 'cigncf', 'cmhypfl', 
           'nadirgp2', 'cd4', 'fgluc')

################################################################################
# PREP DATA FOR CASUAL SURVIVAL FOREST 
################################################################################

# Create dataframe of all covariates in original form
# All covariates are discrete
Z <- data.frame(adgrfcc[,covar])
head(Z)

# Convert all covariates to factors
factors <- lapply(Z, FUN=as.factor)

# Convert categorical variables to one-hot encoded/dummy variables
# fullRank=F keeps all dummy variables => can later choose reference level
dummy <- dummyVars(~ ., data=factors, fullRank=FALSE)

# Construct dummy variable data frame for modeling
X <- data.frame(predict(dummy, newdata=Z))

# Remove reference levels
ref <- c('natalsex.Female', 'gbdrace.HighIncomeWhite', 'cigncf.Never', 'cmhypfl.N',
         'nadirgp2.GTE200')

X <- X %>% select(-all_of(ref))

# Treatment indicator (W=0: Pitavastatin, W=1: Placebo)
# Reverse coded so treatment effects in expected direction (negative = risk red with statin)
W <- ifelse(adgrfcc$rxl == 'Pitavastatin', 0, 1) %>% as.vector()

# Outcome (Time-to-Confirmed Diabetes, years)
Y <- adgrfcc$avald %>% as.vector()

# Event indicator
D <- adgrfcc$evntd %>% as.vector()

# Combine into one dataframe for plotting, subsets, etc
DF <- data.frame(X=X, Y=Y, W=W, D=D)

# Number of rows (participants), Number of columns (variables)
dim(DF)

################################################################################
# FIT CASUAL SURVIVAL FOREST + EXPORT
################################################################################

# number of variables
p = ncol(X)

# Tunable parameters
crf.num.trees     = 10000         # the more trees the better but longer runtime
crf.min.node.size = 30            # Small node size suggested due to rare events
crf.alpha         = 0.01         # Small alpha suggested due to rare events

# Fit causal survival forest
# NOTE 1: Warnings can be ignored. Estimates suggest no evidence of instability,
# and computational time is not a concern. Six year horizon provides best trade-off,
# even if some censoring weights are small. 
# NOTE 2: Censoring weight instability was not an issue for categorical risk factor CSF,
# likely due to fewer strata within which to estimate weights.
crf <- causal_survival_forest(X, Y, W, D, 
                              # (Restricted) survival parameter
                              # ATE = difference in survival probability at year 6
                              target='survival.probability',
                              # time horizon
                              # yr 6 = good tradeoff b/w events and stable IPCW
                              horizon=5,
                              # Treatment propensity (known, 1:1 randomization)
                              W.hat=0.5,
                              # honest splitting for valid inference
                              # keeping in mind rare events
                              honesty=TRUE,
                              honesty.fraction = 0.5,
                              # Number of trees for forest
                              num.trees=crf.num.trees,
                              # small leaves risk overfit, large risk underfit
                              alpha=crf.alpha,
                              sample.fraction=0.5,
                              min.node.size=crf.min.node.size,
                              seed=123)

# Save fit CRF so it can be used for later analyses (SLOW)
saveRDS(crf, file=paste0(rdata, "crf_diab.RDS"))

################################################################################
# PREDICT INDIVIDUAL TREATMENT EFFECTS (ITEs) + EXPORT
################################################################################

# Set estimate.variance=T to be able to construct CIs
pred <- predict(crf, estimate.variance=T)

# Add patids for later joins
ite <- data.frame(patid    = adgrfcc$patid, 
                  ite_diab = pred$predictions, 
                  se_diab  = sqrt(pred$variance.estimates))

# Add patid to covariate design matrix
X$patid <- adgrfcc$patid

# Join ITE predictions and covariate design matrix for later reference
OUT <- ite %>% left_join(X, by='patid')

# Add labels
attr(OUT,          "label") <- "Predicted Individual Treatment Effects and Design Matrix: Diabetes"
attr(OUT,        "comment") <- paste("Created by:", this.path())
attr(OUT$patid,    'label') <- "Participant ID"
attr(OUT$ite_diab, 'label') <- "Predicted ITE: 5-yr Diabetes"
attr(OUT$se_diab,  'label') <- "SE of ITE: 5-yr Diabetes"

# Save to RDS
saveRDS(OUT, file=paste0(rdata, 'ite_diab.RDS'))

# Create "contents" file
write_contents(df=OUT, file=paste0(rdata, 'contents.ite_diab.txt'))
