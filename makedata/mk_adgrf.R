################################################################################
# PROD PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/makedata/production/mk_adgrf.R
# WORK PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/makedata/working/mk_adgrf.R
# 
# RISK ASSESSMENT: MODERATE
# DOCUMENT ID:     PROG.10084.v3
# 
# PURPOSE:         Convert adgrf.sas7bdat to adgrf.RDS for use in R
# 
# SOURCE PRGM:     none
# INPUT:           reprieve/POOLED/final/MACE_hte/sas_datasets/derived/adgrf.sas7bdat
# OUTPUT:          reprieve/POOLED/final/MACE_hte/sas_datasets/Rderived/adgrf.RDS
#
# PACKAGES USED:   this.path, tidyverse, skimr
# FUNCTIONS USED:  read.sas, write_contents
# EXEMPTIONS:      none
# 
# AUTHOR:          Emma Davies Smith
# CREATION DATE:   09FEB2026
# 
# NOTES:           Also adds flags for complete cases (needed for modelling)
#                  Removes participants with pre-existing diabetes (diabfl=Y)
# MODIFICATIONS:   [eds 06apr2026] added write_contents from RWG
################################################################################

################################################################################
# SETUP
################################################################################

# Project setupfile
source("/home/reprieve/POOLED/final/MACE_hte/programs/setupfile.R")

# CBAR's read.sas function
source("/home/programs/R/read.sas.R")

# Library to extract paths
library(this.path ,lib=pkgdir)

# RWG's function to mimick contents.lst (requires tidyverse and skimr)
source("/home/reprieve/POOLED/final/MACE_hte/functions/write_contents.R")
library(tidyverse)
library(skimr, lib=pkgdir)

################################################################################
# LOAD
################################################################################

# Import adgrf SAS dataset
adgrf <- read.sas(libname=sdata, dsname="adgrf")

# View dataframe attributes
str(adgrf)

# Preview adgrf
head(adgrf, 20)

# Exclude participants with pre-existing diabetes
adgrf <- adgrf[adgrf$diabfl!='Y',]

# Check missing data
apply(adgrf, MARGIN=2, FUN=function(x) sum(is.na(x)))

################################################################################
# EDIT
################################################################################

# Add flag for complete cases
covar <- c('gbdrace', 'natalsex', 'age', 
           'sysbp', 'schol', 'shdl', 'bmi', 'egfr_eq21', 'cigncf', 'cmhypfl', 
           'nadirgp2', 'cd4', 'fgluc')

adgrf$ccfl <- complete.cases(adgrf[,c('patid', covar)])

################################################################################
# LABEL
################################################################################

# Add variable labels and units 
attr(adgrf$patid     ,'label')  <- "Participant ID"
attr(adgrf$rxl       ,'label')  <- "Randomized arm"
attr(adgrf$rxln      ,'label')  <- "Randomized arm (numeric)"
attr(adgrf$pce       ,'label')  <- "10-year PCE ASCVD risk score"
attr(adgrf$pce       ,'units')  <- "%"
attr(adgrf$pcefl     ,'label')  <- "PCE flag (<5%)"
attr(adgrf$bool3     ,'label')   <- "Low risk boolean"
attr(adgrf$age       ,'label')  <- "Age" 
attr(adgrf$age       ,'units')  <- "years"
attr(adgrf$agegp     ,'label')  <- "Age"
attr(adgrf$agegp     ,'units')  <- "category"
attr(adgrf$agegp3    ,'label')  <- "Age"
attr(adgrf$natalsex  ,'label')  <- "Natal sex"
attr(adgrf$racegp    ,'label')  <- "Race"
attr(adgrf$cntryg1   ,'label')  <- "Country"
attr(adgrf$gbdgp1    ,'label')  <- "GBD Region"
attr(adgrf$gbdrace   ,'label')  <- "GBD region, Race in HI"
attr(adgrf$cmdiabfl  ,'label')  <- "On anti-diabetic medication"
attr(adgrf$cmhypfl   ,'label')  <- "On anti-hypertensive medication"
attr(adgrf$schol     ,'label')  <- "Total cholesterol"
attr(adgrf$diabp     ,'label')  <- "Diastolic blood pressure"
attr(adgrf$diabp     ,'units')  <- "mmHg"
attr(adgrf$egfr_eq21 ,'label')  <- "eGFR"
attr(adgrf$egfr_eq21 ,'units')  <- "mL/min/1.73m\u00B2"
attr(adgrf$egfrfl    ,'label')  <- "eGFR"
attr(adgrf$egfrfl    ,'units')  <- "category"
attr(adgrf$cigstgp   ,'label')  <- "Smoking status"
attr(adgrf$subsgp1   ,'label')  <- "Substance use"
attr(adgrf$bmi       ,'label')  <- "BMI"
attr(adgrf$bmigp6    ,'units')  <- "kg/m\u00B2"
attr(adgrf$sysbp     ,'label')  <- "Systolic blood pressure"
attr(adgrf$sysbp     ,'units')  <- "mmHg"
attr(adgrf$shdl      ,'label')  <- "HDL-C"
attr(adgrf$shdl      ,'units')  <- "mg/dL"
attr(adgrf$fgluc     ,'label')  <- "Fasting glucose"
attr(adgrf$fgluc     ,'units')  <- "mg/dL"  
attr(adgrf$cnsrm     ,'label')  <- "Censor indicator: MACE"
attr(adgrf$evntm     ,'label')  <- "Event indicator: MACE"
attr(adgrf$avalm     ,'label')  <- "Time-to-event: MACE"
attr(adgrf$cnsrd     ,'label')  <- "Censor indicator: Diabetes"
attr(adgrf$evntd     ,'label')  <- "Event indicator: Diabetes"
attr(adgrf$avald     ,'label')  <- "Time-to-event: Diabetes"
attr(adgrf$ccfl      ,'label')  <- "Complete case flag"
attr(adgrf$nadirgp2  ,'label')  <- "Nadir CD4"
attr(adgrf$nadirgp2  ,'units')  <- "cells/mm\u00B3"
attr(adgrf$cd4       ,'label')  <- "CD4 count"
attr(adgrf$cd4       ,'units')  <- "cells/mm\u00B3"
attr(adgrf$pvt_ascvd10,'label')   <- "PREVENT-ASCVD risk score"
attr(adgrf$pvt_ascvd10g1,'label') <- "PREVENT-ASCVD risk score, categories"
attr(adgrf$pvt_lowfl,'label')     <- "Risk strata"

# Set path label for traceability
attr(adgrf, "label")    <- "Generalized Random Forest (GRF) Inputs"
attr(adgrf, "comment")  <- paste("Created by:", this.path())
################################################################################
# SAVE
################################################################################

# Save to RDS
saveRDS(adgrf, file=paste0(rdata, 'adgrf.RDS'))

# Use write_contents script written by RWG
# Requires skimr (and tidyverse in modified version)
write_contents(df=adgrf, file=paste0(rdata, 'contents.adgrf.txt'))
               
               