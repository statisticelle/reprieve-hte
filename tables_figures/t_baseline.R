################################################################################
# PROD PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/production/t_baseline.R
# WORK PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/working/t_baseline.R
# 
# RISK ASSESSMENT: Low
# DOCUMENT ID:     PROG.10084.v3
# 
# PURPOSE:         Summarize baseline characteristics of REPRIEVE participants
#                  that are included in survival forest and analysis
# 
# SOURCE PRGM:     none
# INPUT:           /home/reprieve/POOLED/final/MACE_hte/sas_datasets/Rderived/adgrf.RDS
# OUTPUT:          /home/reprieve/POOLED/final/MACE_hte/output/tables/t_baseline1_mace.rtf
#                                                                     t_baseline2_mace.rtf
#                                                                     t_baseline1_diab.rtf
#                                                                     t_baseline2_diab.rtf
#
# PACKAGES USED:   tidyverse, flextable, officer, table1, this.path
# FUNCTIONS USED:  r_word_helpers.R, rtf_numbering.R          
# EXEMPTIONS:      none
# 
# AUTHOR:          Emma Davies Smith
# CREATION DATE:   19NOV2025
# 
# NOTES:           Only a subset of covariates are reported (Those used in modelling)
# MODIFICATIONS:   [eds 07apr2026] removed croi; added cts vars to baseline tables
################################################################################

################################################################################
# SETUP
################################################################################

# Project setupfile
source("/home/reprieve/POOLED/final/MACE_hte/programs/setupfile.R")

# Libraries
# NOTE: These packages will not install within the project dir due to issue with systemfonts
#       I was told to use RWG versions as work-around (lib=rwgdir)
library(tidyverse) 
library(flextable ,lib=rwgdir)
library(officer   ,lib=rwgdir)
library(table1    ,lib=pkgdir)
library(this.path ,lib=pkgdir)

# Reporting helper functions (from RWG)
source("/home/reprieve/POOLED/final/MACE_hte/functions/r_word_helpers.R")
source("/home/reprieve/POOLED/final/MACE_hte/functions/rtf_numbering.R")

# Setup page size
pr <- prop_section(page_size = page_size(orient='portrait', width=8.5, height=11),
                   page_margins = page_mar(bottom=0.5, top=0.5, left=0.5, right=0.5, 
                                           header=0.25, footer=0.25))

# Set working directory
setwd(this.dir())

################################################################################
# IMPORT DATA
################################################################################

# Load dataset: adgrf
adgrf <- readRDS(paste0(rdata, 'adgrf.RDS'))

# Restrict to complete cases
adgrf <- adgrf %>% filter(ccfl)

################################################################################
# Event status
################################################################################

adgrf$evntm5 <- ifelse(adgrf$evntm == 1 & adgrf$avalm < 5, 'MACE', 
                 ifelse(adgrf$evntm == 0 & adgrf$avalm >= 5, 'No MACE',
                        'Censored'))
adgrf$evntd5 <- ifelse(adgrf$evntd == 1 & adgrf$avald < 5, 'Diabetes', 
                 ifelse(adgrf$evntd == 0 & adgrf$avald >= 5, 'No diabetes',
                      'Censored'))

################################################################################
# RELABEL AND REORDER FACTORS
################################################################################

# Format/order categorical values
adgrf$gbdrace <- factor(adgrf$gbdrace, 
                        levels=c("HighIncomeWhite", "HighIncomeBlack", "HighIncomeOther", 
                                 "LatinAmericaandCaribbean", "S.East/EastAsia",
                                 "SouthAsia", "Sub-SaharanAfrica"),
                        labels=c("High Income, White", "High Income, Black", "High Income, Other",
                                 "Latin America and Caribbean", "S.East/East Asia", 
                                 "South Asia", "Sub-Saharan Africa"))

adgrf$cmhypfl   <- factor(adgrf$cmhypfl, levels=c('Y', 'N'), labels=c('Yes', 'No'))
adgrf$nadirgp2  <- factor(adgrf$nadirgp2, levels=c('LT200', 'GTE200'), labels=c('<200', '\u2265 200'))
adgrf$evntm5    <- factor(adgrf$evntm5, levels=c('MACE',     'No MACE', 'Censored'))
adgrf$evntd5    <- factor(adgrf$evntd5, levels=c('Diabetes', 'No diabetes', 'Censored'))

adgrf$pvt_ascvd10g1 <- factor(adgrf$pvt_ascvd10g1, levels=c('0-<3', '3-<5', '5-<10', '10+'),
                              labels=c('0-<3', '3-<5', '5+', '5+'))
################################################################################
# REAPPLY LABELS
################################################################################

# Add variable labels and units 
attr(adgrf$patid     ,'label')  <- "Participant ID"
attr(adgrf$rxl       ,'label')  <- "Randomized arm"
attr(adgrf$rxln      ,'label')  <- "Randomized arm (numeric)"
attr(adgrf$age       ,'label')  <- "Age" 
attr(adgrf$age       ,'units')  <- "years"
attr(adgrf$agegp     ,'label')  <- "Age"
attr(adgrf$agegp     ,'units')  <- "years"
attr(adgrf$agegp3    ,'label')  <- "Age"
attr(adgrf$natalsex  ,'label')  <- "Natal sex"
attr(adgrf$gbdrace   ,'label')  <- "GBD region, Race"
attr(adgrf$cmhypfl   ,'label')  <- "Anti-hypertensive medication"
attr(adgrf$egfr_eq21 ,'label')  <- "eGFR"
attr(adgrf$egfr_eq21 ,'units')  <- "mL/min/1.73m\u00B2"
attr(adgrf$cigncf    ,'label')  <- "Smoking status"
attr(adgrf$bmi       ,'label')  <- "BMI"
attr(adgrf$bmi       ,'units')  <- "kg/m\u00B2"
attr(adgrf$sysbp     ,'label')  <- "Systolic blood pressure"
attr(adgrf$sysbp     ,'units')  <- "mmHg"
attr(adgrf$shdl      ,'label')  <- "HDL cholesterol"
attr(adgrf$shdl      ,'units')  <- "mg/dL"
attr(adgrf$schol     ,'label')  <- "Total cholesterol"
attr(adgrf$schol     ,'units')  <- "mg/dL"
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
attr(adgrf$pvt_ascvd10,  'label') <- "PREVENT-ASCVD score"
attr(adgrf$pvt_ascvd10,  'units') <- "%"
attr(adgrf$pvt_ascvd10g1,'label') <- "PREVENT-ASCVD score"
attr(adgrf$pvt_ascvd10g1,'units') <- "%"
attr(adgrf$evntm5,       'label') <- "Event status by 5 years of follow-up"
attr(adgrf$avalm,        'label') <- "Follow-up"
attr(adgrf$avalm,        'units') <- "years"
attr(adgrf$avald,        'label') <- "Follow-up"
attr(adgrf$avald,        'units') <- "years"

################################################################################
# BASELINE DEMOGRAPHICS AND EXPOSURES: MACE
################################################################################

# Create a caption object for use in a flextable
# * turn this into a helper function
cap1 <- as_paragraph(as_chunk('Randomized arm, sociodemographics, and risk scores by MACE event status', 
                             props=fp_text(font.size=my.caption.font.size, font.family='Arial', bold=T)))
          
# Set numbering and bookmark for table
# * this is a much simpler version than whats used by sean's demogs_tbl.R
auto1 <- run_autonum(seq_id="Table", bkm="baseline1_mace", 
                    prop=fp_text(font.size=my.caption.font.size, font.family='Arial', bold=T))

# Create table1 object
# NOTE: This will yield warning if missing values for only one arm (ignorable)
tab1 <- table1(
  x = ~ rxl + gbdrace + age + natalsex + pvt_ascvd10 + pvt_ascvd10g1 | evntm5,
  data = adgrf,
  render.categorical = render.categorical.nomiss,
  render.continuous  = render.continuous.nomiss.detailed,
  render.missing     = NULL
) %>% suppressWarnings()

# Conver table1 object to flextable and format
ftab1 <- tab1 %>% 
          t1flex() %>% 
          paginate_by_variable() %>%
          border_inner_h(border=fp_border(col="gray92")) %>%
          height(part='body', height=0.25) %>%
          hrule(rule='exact', part='body') %>%
          height(part='header', height=0.45) %>%
          hrule(rule='exact', part='header') %>%
          width(j=1, width=2.25) %>%
          width(j=2:5, width=1) %>%
          padding(part='all', padding=0) %>%
          line_spacing(part='all', space=1) %>%
          font(part='all', fontname='Arial') %>%
          fontsize(part='all', size=8) %>%
          valign(part='all', valign='center') %>%
          pgmfootdt(path=this.path(), size=my.foot.font.size) %>%
          set_caption(caption=cap1, autonum=auto1) %>%
          set_table_properties(layout='fixed')

ftab1

save_as_rtf(ftab1, path=paste0(rtab, 't_baseline1_mace.rtf'), pr_section=pr)
update_rtf_table_numbering(paste0(rtab, 't_baseline1_mace.rtf'), overwrite=T)

################################################################################
# BASELINE CONTINUED
################################################################################

# Create a caption object for use in a flextable
# * turn this into a helper function
cap2 <- as_paragraph(as_chunk('Baseline risk factors by MACE event status', 
                              props=fp_text(font.size=my.caption.font.size, font.family='Arial', bold=T)))

# Set numbering and bookmark for table
# * this is a much simpler version than whats used by sean's demogs_tbl.R
auto2 <- run_autonum(seq_id="Table", bkm="baseline2_mace", 
                     prop=fp_text(font.size=my.caption.font.size, font.family='Arial', bold=T))

# Create table1 object
# NOTE: This will yield warning if missing values for only one arm (ignorable)
tab2 <- table1(
  x = ~ sysbp + schol + shdl + bmi + egfr_eq21 + cigncf + cmhypfl + nadirgp2 + cd4 + fgluc | evntm5,
  data = adgrf,
  render.categorical = render.categorical.nomiss,
  render.continuous  = render.continuous.nomiss.detailed,
  render.missing     = NULL
) %>% suppressWarnings()

# Conver table1 object to flextable and format
ftab2 <- tab2 %>% 
  t1flex() %>% 
  paginate_by_variable() %>%
  border_inner_h(border=fp_border(col="gray92")) %>%
  height(part='body', height=0.25) %>%
  hrule(rule='exact', part='body') %>%
  height(part='header', height=0.45) %>%
  hrule(rule='exact', part='header') %>%
  width(j=1, width=2.25) %>%
  width(j=2:5, width=1) %>%
  padding(part='all', padding=0) %>%
  line_spacing(part='all', space=1) %>%
  font(part='all', fontname='Arial') %>%
  fontsize(part='all', size=8) %>%
  valign(part='all', valign='center') %>%
  pgmfootdt(path=this.path(), size=my.foot.font.size) %>%
  set_caption(caption=cap2, autonum=auto2) %>%
  set_table_properties(layout='fixed')

ftab2

save_as_rtf(ftab2, path=paste0(rtab, 't_baseline2_mace.rtf'), pr_section=pr)
update_rtf_table_numbering(paste0(rtab, 't_baseline2_mace.rtf'), overwrite=T)

################################################################################
# BASELINE DEMOGRAPHICS AND EXPOSURES: DIABETES
################################################################################

# Create a caption object for use in a flextable
# * turn this into a helper function
cap1 <- as_paragraph(as_chunk('Randomized arm, sociodemographics, and risk scores by DM event status', 
                              props=fp_text(font.size=my.caption.font.size, font.family='Arial', bold=T)))

# Set numbering and bookmark for table
# * this is a much simpler version than whats used by sean's demogs_tbl.R
auto1 <- run_autonum(seq_id="Table", bkm="baseline1_diab", 
                     prop=fp_text(font.size=my.caption.font.size, font.family='Arial', bold=T))

# Create table1 object
# NOTE: This will yield warning if missing values for only one arm (ignorable)
tab1 <- table1(
  x = ~ rxl + gbdrace + age + natalsex + pvt_ascvd10 + pvt_ascvd10g1 | evntd5,
  data = adgrf,
  render.categorical = render.categorical.nomiss,
  render.continuous  = render.continuous.nomiss.detailed,
  render.missing     = NULL
) %>% suppressWarnings()

# Conver table1 object to flextable and format
ftab1 <- tab1 %>% 
  t1flex() %>% 
  paginate_by_variable() %>%
  border_inner_h(border=fp_border(col="gray92")) %>%
  height(part='body', height=0.25) %>%
  hrule(rule='exact', part='body') %>%
  height(part='header', height=0.45) %>%
  hrule(rule='exact', part='header') %>%
  width(j=1, width=2.25) %>%
  width(j=2:5, width=1) %>%
  padding(part='all', padding=0) %>%
  line_spacing(part='all', space=1) %>%
  font(part='all', fontname='Arial') %>%
  fontsize(part='all', size=8) %>%
  valign(part='all', valign='center') %>%
  pgmfootdt(path=this.path(), size=my.foot.font.size) %>%
  set_caption(caption=cap1, autonum=auto1) %>%
  set_table_properties(layout='fixed')

ftab1

save_as_rtf(ftab1, path=paste0(rtab, 't_baseline1_diab.rtf'), pr_section=pr)
update_rtf_table_numbering(paste0(rtab, 't_baseline1_diab.rtf'), overwrite=T)

################################################################################
# BASELINE CONTINUED
################################################################################

# Create a caption object for use in a flextable
# * turn this into a helper function
cap2 <- as_paragraph(as_chunk('Baseline risk factors by DM event status', 
                              props=fp_text(font.size=my.caption.font.size, font.family='Arial', bold=T)))

# Set numbering and bookmark for table
# * this is a much simpler version than whats used by sean's demogs_tbl.R
auto2 <- run_autonum(seq_id="Table", bkm="baseline2_diab", 
                     prop=fp_text(font.size=my.caption.font.size, font.family='Arial', bold=T))

# Create table1 object
# NOTE: This will yield warning if missing values for only one arm (ignorable)
tab2 <- table1(
  x = ~ sysbp + schol + shdl + bmi + egfr_eq21 + cigncf + cmhypfl + nadirgp2 + cd4 + fgluc | evntd5,
  data = adgrf,
  render.categorical = render.categorical.nomiss,
  render.continuous  = render.continuous.nomiss.detailed,
  render.missing     = NULL
) %>% suppressWarnings()

# Conver table1 object to flextable and format
ftab2 <- tab2 %>% 
  t1flex() %>% 
  paginate_by_variable() %>%
  border_inner_h(border=fp_border(col="gray92")) %>%
  height(part='body', height=0.25) %>%
  hrule(rule='exact', part='body') %>%
  height(part='header', height=0.45) %>%
  hrule(rule='exact', part='header') %>%
  width(j=1, width=2.25) %>%
  width(j=2:5, width=1) %>%
  padding(part='all', padding=0) %>%
  line_spacing(part='all', space=1) %>%
  font(part='all', fontname='Arial') %>%
  fontsize(part='all', size=8) %>%
  valign(part='all', valign='center') %>%
  pgmfootdt(path=this.path(), size=my.foot.font.size) %>%
  set_caption(caption=cap2, autonum=auto2) %>%
  set_table_properties(layout='fixed')

ftab2

save_as_rtf(ftab2, path=paste0(rtab, 't_baseline2_diab.rtf'), pr_section=pr)
update_rtf_table_numbering(paste0(rtab, 't_baseline2_diab.rtf'), overwrite=T)

