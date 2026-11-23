################################################################################
# PROD PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/production/f_partialdep.r
# WORK PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/working/f_partialdep.r
# 
# RISK ASSESSMENT: Low
# DOCUMENT ID:     PROG.10084.v3
# 
# PURPOSE:         Partial dependence plots
# 
# SOURCE PRGM:     none
# INPUT:           
# OUTPUT:          
#
# PACKAGES USED:   none
# FUNCTIONS USED:  none
# EXEMPTIONS:      none
# 
# AUTHOR:          Emma Davies Smith
# CREATION DATE:   04DEC2025
# 
# NOTES:           none
# MODIFICATIONS:   none 
################################################################################


################################################################################
# SETUP
################################################################################

# Project setupfile
source("/home/reprieve/POOLED/final/MACE_hte/programs/setupfile.R")

# Analysis libraries
library(tidyverse)              # data manipulation
library(grf,       lib=pkgdir)  # generalized random forests
library(lmtest,    lib=pkgdir)  # dependency of grf that doesnt load properly
library(sandwich,  lib=pkgdir)  # dependency of grf that doesnt load properly 
library(survival)               # survival methods (default package)
library(caret,     lib=pkgdir)  # create dummy variables
library(gridExtra, lib=pkgdir)

# Reporting libraries
library(flextable ,lib=rwgdir)
library(officer   ,lib=rwgdir)
library(table1    ,lib=rwgdir)
library(this.path ,lib=rwgdir)

# Reporting helper functions (from RWG)
source("/home/reprieve/POOLED/final/MACE_hte/functions/r_word_helpers.R")
source("/home/reprieve/POOLED/final/MACE_hte/functions/rtf_numbering.R")

################################################################################
# LOAD DATA 
################################################################################

# Load dataset: adgrf
adgrf <- readRDS(paste0(rdata, 'adgrf.RDS'))

# Subset to complete cases
adgrfcc <- adgrf %>% filter(ccfl)

# Load predicted ITEs
ite_mace <- readRDS(paste0(rdata, 'ite_mace.RDS'))
ite_diab <- readRDS(paste0(rdata, 'ite_diab.RDS'))

ite <- adgrfcc %>%
  left_join(ite_mace %>% select(patid, ite_mace)) %>%
  left_join(ite_diab %>% select(patid, ite_diab))

################################################################################
# LABEL DATA
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
attr(adgrf$bmigp6    ,'label')  <- "Race-specific BMI group"
attr(adgrf$bmigp6    ,'units')  <- "kg/m\u00B2"
attr(adgrf$diabfl    ,'label')  <- "Pre-existing diabetes"
attr(adgrf$hypfl     ,'label')  <- "Pre-existing hypertension"
attr(adgrf$wstcir    ,'label')  <- "Waist circumfrence"
attr(adgrf$wstcir    ,'units')  <- "cm"
attr(adgrf$ms_wcir   ,'label')  <- "Waist circumfrence"
attr(adgrf$ms_wcir   ,'units')  <- "category"
attr(adgrf$metS      ,'label')  <- "Metabolic syndrome"
attr(adgrf$metS2     ,'label')  <- "Metabolic syndrome (indeterminate=missing)"
attr(adgrf$fhcvd     ,'label')  <- "Family history of CVD"
attr(adgrf$fhcvd2    ,'label')  <- "Family history of CVD (unknown=missing)"
attr(adgrf$prediab   ,'label')  <- "Prediabetic status"
attr(adgrf$prediab2  ,'label')  <- "Prediabetic status (indeterminate=missing)"
attr(adgrf$sysbp     ,'label')  <- "Systolic blood pressure"
attr(adgrf$sysbp     ,'units')  <- "mmHg"
attr(adgrf$hscrp     ,'label')  <- "hs-CRP"
attr(adgrf$hscrp     ,'units')  <- "mg/L"
attr(adgrf$shdl      ,'label')  <- "HDL-C"
attr(adgrf$shdl      ,'units')  <- "mg/dL"
attr(adgrf$hdlfl     ,'label')  <- "HDL-C"
attr(adgrf$hdlfl     ,'units')  <- "category"
attr(adgrf$sldlc     ,'label')  <- "LDL-C"
attr(adgrf$sldlc     ,'units')  <- "mg/dL"
attr(adgrf$ldlfl     ,'label')  <- "LDL-C"
attr(adgrf$ldlfl     ,'units')  <- "category"
attr(adgrf$strig     ,'label')  <- "Triglycerides"
attr(adgrf$strig     ,'units')  <- "mg/dL"  
attr(adgrf$trigfl    ,'label')  <- "Triglycerides"
attr(adgrf$trigfl    ,'units')  <- "category"
attr(adgrf$fgluc     ,'label')  <- "Fasting glucose"
attr(adgrf$fgluc     ,'units')  <- "mg/dL"  
attr(adgrf$fglucfl   ,'label')  <- "Fasting glucose"
attr(adgrf$fglucfl   ,'units')  <- "category"
attr(adgrf$cnsrm     ,'label')  <- "Censor indicator: MACE"
attr(adgrf$evntm     ,'label')  <- "Event indicator: MACE"
attr(adgrf$avalm     ,'label')  <- "Time-to-event: MACE"
attr(adgrf$cnsrd     ,'label')  <- "Censor indicator: Diabetes"
attr(adgrf$evntd     ,'label')  <- "Event indicator: Diabetes"
attr(adgrf$avald     ,'label')  <- "Time-to-event: Diabetes"
attr(adgrf$ccfl      ,'label')  <- "Complete case flag"
attr(adgrf$nadirgp2  ,'label')  <- "Nadir CD4"
attr(adgrf$nadirgp2  ,'units')  <- "category"
attr(adgrf$cd4       ,'label')  <- "CD4 count"
attr(adgrf$cd4gp3    ,'label')  <- "CD4 count group"
attr(adgrf$cd4gp3    ,'units')  <- "category"
attr(adgrf$arvdur    ,'label')  <- "Lifetime ART duration"
attr(adgrf$arvdurgp  ,'label')  <- "Lifetime ART duration"
attr(adgrf$arvdurgp  ,'units')  <- "category"
attr(adgrf$acat      ,'label')  <- "Baseline ART regimen"
attr(adgrf$acat2     ,'label')  <- "Baseline ART regimen"
attr(adgrf$abcfl     ,'label')  <- "Current or prior abacavir use"
attr(adgrf$instifl   ,'label')  <- "Current INSTI use"
attr(adgrf$hiv1rnaq  ,'label')  <- "Quantifiable HIV-1 RNA"
attr(adgrf$hiv1rnagp2,'label')  <- "HIV-1 RNA group"
attr(adgrf$pvt_ascvd10,'label')   <- "PREVENT-ASCVD risk score"
attr(adgrf$pvt_ascvd10g1,'label') <- "PREVENT-ASCVD risk score, categories"
attr(adgrf$pvt_lowfl,'label')     <- "Risk strata"

ite$racegp  <- factor(ite$racegp,  
                      levels=c('Black or African American', 'White', 'Asian', 'Other'))
ite$gbdrace <- factor(ite$gbdrace, 
                      levels=c("HighIncomeWhite", "HighIncomeBlack", "HighIncomeOther",
                               "LatinAmericaandCaribbean", "S.East/EastAsia",
                               "SouthAsia", "Sub-SaharanAfrica"),
                      labels=c("HI-W", "HI-B", "HI-O",
                               "LatAmCa", "SE/EAsia", 
                               "SAsia", "SSAfrica"))
ite$natalsex <- factor(ite$natalsex, levels=c('Female', 'Male'))
ite$cigncf   <- factor(ite$cigncf, levels=c('Never', 'Former', 'Current'))
ite$nadirgp2  <- factor(ite$nadirgp2, levels=c('LT200', 'GTE200'), labels=c('<200', '\u2265 200'))
ite$cmhypfl <- factor(ite$cmhypfl, levels=c('N', 'Y'), labels=c('No', 'Yes'))

# Need to clean up variable labels for reporting
# This maps the label from BLP to variable names used by baseline tables
label_map <- c("(Intercept)" = "(Intercept)",
               "age"         ="Age (years)",
               "natalsex"    ="Natal sex",
               "gbdrace"     ="GBD region, Race",
               "cigncf"      ="Smoking status",
               "cmhypfl"     ="Anti-hypertensive medication",
               "bmi"         ="BMI (kg/m\u00B2)",
               "sysbp"       ="Systolic blood pressure (mmHg)",
               "shdl"        ="HDL cholesterol (mg/dL)",
               "schol"       ="Total cholesterol (mg/dL)",
               "trigfl"      ="Triglycerides (mg/dL)",
               "fgluc"      ="Fasting glucose (mg/dL)",
               "egfr_eq21"   ="eGFR (mL/min/1.73m\u00B2)",
               "nadirgp2"    ="Nadir CD4 (cells/mm\u00B3)",
               "cd4"         ="CD4 count (cells/mm\u00B3)",
               "hscrp"       ="hs-CRP",
               "arvdur"      ="Lifetime ART duration (years)",
               "abcfl"       ="Current or prior abacavir use",
               "instifl"     ="Current INSTI use"
)

# Convert label_map to a dataframe to permit join with formatted data
label_df <- data.frame(short=names(label_map), label=label_map)

################################################################################
# MACE: Partial dependence plots
################################################################################

# Specify covariates of interest
# !!! NOTE: MAKE SURE THIS ALIGNS WITH covar SPECIFIED IN mk_adgrf.R !!!
# NOTE: Treatment arm is not included as part of covar.
covar <- c('gbdrace', 'natalsex', 'age', 
           'sysbp', 'schol', 'shdl', 'bmi', 'egfr_eq21', 'cigncf', 'cmhypfl', 
           'nadirgp2', 'cd4', 'fgluc')

pdp_plot <- function(v, ite_var){
  ymin <- min(ite[[ite_var]])
  ymax <- max(ite[[ite_var]])
  
  if (!is.factor(ite[[v]])) {
    xmin <- quantile(ite[[v]], prob=0.01)
    xmax <- quantile(ite[[v]], prob=0.99)
    ite %>%
      ggplot(aes(x=.data[[v]], y=.data[[ite_var]])) +
        geom_point(alpha=0.15) +
        geom_smooth() +
        theme_bw() +
        theme(axis.text.x=element_text(size=6),
              plot.title=element_text(size=10),
              plot.margin=margin(t=1,r=2,b=1,l=0)) +
        xlim(xmin, xmax) +
        scale_y_continuous(labels=function(x) paste0(x*100, '%'),
                         limits=c(ymin, ymax), 
                         breaks=seq(round(ymin, 2), round(ymax, 2), 0.01)) +
        labs(x=NULL, title=label_map[{{v}}], y='')
  } else if(is.factor(ite[[v]])){
    ite %>%
      ggplot(aes(x=.data[[v]], y=.data[[ite_var]], fill=.data[[v]])) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x=element_text(size=6),
            plot.title=element_text(size=10),
            plot.margin=margin(t=1,r=2,b=1,l=0)) +
      scale_y_continuous(labels=function(x) paste0(x*100, '%'),
                         limits=c(ymin, ymax), 
                         breaks=seq(round(ymin, 2), round(ymax, 2), 0.01)) +
      guides(fill='none') +
      labs(x=NULL, title=label_map[{{v}}], y='')
  } else {
    print('Error')
  }
}

mace_plots <- list()
for(i in 1:length(covar)){
  mace_plots[[i]] <- pdp_plot(v=covar[i], ite_var='ite_mace')
}

grid_mace <- grid.arrange(grobs=mace_plots, ncol=2)

ggsave(paste0(rfig, 'f_partialdep/pdp_mace.png'), grid_mace, width=8.5, height=11)

insertPicture_rtf(plotin=paste0(rfig, 'f_partialdep/pdp_mace.png'), 
                  path=this.path(),
                  fileout_loc=rfig, 
                  caption.font.size=my.caption.font.size,
                  foot.font.size=my.foot.font.size,
                  title="Predicted Conditional Average Treatment Effects (CATEs) for MACE by risk factor", 
                  bkm="pdp_mace", 
                  width=8.5*0.7, height=11*0.7)

update_rtf_fig_numbering(paste0(rfig, 'pdp_mace.rtf'), overwrite=T)

diab_plots <- list()
for(i in 1:length(covar)){
  diab_plots[[i]] <- pdp_plot(v=covar[i], ite_var='ite_diab')
}

grid_diab <- grid.arrange(grobs=diab_plots, ncol=2)

ggsave(paste0(rfig, 'f_partialdep/pdp_diab.png'), grid_diab, width=8.5, height=11)

insertPicture_rtf(plotin=paste0(rfig, 'f_partialdep/pdp_diab.png'), 
                  path=this.path(),
                  fileout_loc=rfig, 
                  caption.font.size=my.caption.font.size,
                  foot.font.size=my.foot.font.size,
                  title="Predicted Conditional Average Treatment Effects (CATEs) for DM by risk factor", 
                  bkm="pdp_diab", 
                  width=8.5*0.7, height=11*0.7)

update_rtf_fig_numbering(paste0(rfig, 'pdp_diab.rtf'), overwrite=T)
