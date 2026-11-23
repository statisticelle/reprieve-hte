################################################################################
# PROD PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/production/t_var_imp.r
# WORK PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/working/t_var_imp.r
# 
# RISK ASSESSMENT: Low
# DOCUMENT ID:     PROG.10084.v3
# 
# PURPOSE:         Causal survival forest variable importance for both events
# 
# SOURCE PRGM:     none
# INPUT:           
#                  rdata/crf_mace.RDS       (causal survival forest object for MACE)
#                  rdata/ite_mace.RDS       (ite predictions and design matrix for MACE)
#                  rdata/crf_diab.RDS       (causal survival forest object for DM)
#                  rdata/ite_diab.RDS       (ite predictions and design matrix for DM)
# OUTPUT:          output/figures/f_var_imp/vimp_mace.png
#                  output/figures/f_var_imp/vimp_diab.png
#                  output/figures/vimp_mace.rtf
#                  output/figures/vimp_diab.rtf
#
# PACKAGES USED:   tidyverse, grf, lmtest, sandwich, survival, caret, flextable, 
#                  officer, table1, this.path
# FUNCTIONS USED:  r_word_helpers.R, rtf_numbering.R
# EXEMPTIONS:      none
# 
# AUTHOR:          Emma Davies Smith
# CREATION DATE:   03DEC2025
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

# Reporting libraries
library(flextable ,lib=rwgdir)
library(officer   ,lib=rwgdir)
library(table1    ,lib=rwgdir)
library(this.path ,lib=rwgdir)

# Reporting helper functions (from RWG)
source("/home/reprieve/POOLED/final/MACE_hte/functions/r_word_helpers.R")
source("/home/reprieve/POOLED/final/MACE_hte/functions/rtf_numbering.R")

################################################################################
# LOAD DATA AND OBJECTS
################################################################################

# Load crf objects (slow)
crf_mace <- readRDS(paste0(rdata, 'crf_mace.RDS'))
crf_diab <- readRDS(paste0(rdata, 'crf_diab.RDS')) 

# Load ITE predictions and design matrices
ite_mace <- readRDS(paste0(rdata, 'ite_mace.RDS'))
ite_diab <- readRDS(paste0(rdata, 'ite_diab.RDS'))

################################################################################
# ASSESS VARIABLE IMPORTANCE: MACE
################################################################################

# Extract variable names used to fit crf
X_mace <- ite_mace %>% select(-all_of(c('patid', 'ite_mace', 'se_mace')))

# Calculate variable importance
vimp_mace <- data.frame(covar=colnames(X_mace),
                        importance=variable_importance(crf_mace)) 

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

# Need to clean up variable labels and values 
# This splits strings, applies regex to cleanup, and formats numeric values
vimp_mace_label <- vimp_mace %>%
  # Variable = first str before .
  # Values   = second, third, fourth str split by .
  mutate(var =str_split_i(covar, '[.]', 1),
         val1=str_split_i(covar, '[.]', 2),
         val2=str_split_i(covar, '[.]', 3),
         val3=str_split_i(covar, '[.]', 4)) %>%
  # Replace empty splits ("NA") with empty str
  mutate_all(~replace_na(., '')) %>%
  # Combine val1, val2, val3 into val and apply regex to get nice value labels
  mutate(val=paste(val1, val2, val3),
         # replace GTE and LT with appropriate symbols
         val=gsub("GTE", "\u2265", val),
         val=gsub("LT", "<", val),
         # add space before capital letter to separate compounded words
         val=trimws(gsub("([A-Z])", " \\1", val)),
         # add hyphen between two numbers
         val=gsub("(\\d) (\\d)", "\\1-\\2", val),
         # add + to standalone number
         val=ifelse(!grepl('-|\u2265|<', val), gsub("(\\d)$", "\\1+", val), val),
         # replace "and" and duplicate spaces with 1 space
         val=gsub("and", " ", val),
         val=gsub("  ", " ", val),
         # replace 'Y' with 'Yes'
         val=ifelse(val=='Y', 'Yes', val),
         val=ifelse(val=='N', 'No', val)) %>%
  left_join(label_df, by=join_by(var==short))  %>%
  # Remove ugly vals and no longer needed covar column
  select(label, var, val, importance) 

vimp_mace_report <- vimp_mace_label %>%
  group_by(label) %>%
  reframe(var=var, val='Total', importance=mean(importance)) %>%
  bind_rows(vimp_mace_label) %>%
  distinct() %>%
  arrange(-importance) %>%
  filter(val=='Total') %>%
  mutate(vid = row_number()) 
  
vimp_mace_report

vimp_mace_plot <- vimp_mace_report %>%
  filter(val=='Total') %>%
  ggplot(aes(y=reorder(label, importance), x=importance)) +
  geom_col(alpha=0.75, aes(fill=label)) +
  geom_line(group=1) +
  geom_point(group=1) +
  labs(x='Variable importance', y=NULL) +
  guides(fill="none") +
  theme_bw() +
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=16, face='bold'),
        axis.title.y=element_text(size=16, face='bold'),
        axis.text.y=element_text(size=16),
        legend.title=element_text(size=12, face='bold'),
        legend.text=element_text(size=12),
        legend.box.background=element_rect(color='black', linewidth=1),
        plot.margin=unit(c(1, 0.5, 5.5, 5.5), 'pt'))

vimp_mace_plot  

# Save to output directory
ggsave(paste0(rfig, 'f_var_imp/vimp_mace.png'), plot=vimp_mace_plot, width=8, height=4.5)

insertPicture_rtf(plotin=paste0(rfig, 'f_var_imp/vimp_mace.png'), 
                  path=this.path(),
                  fileout_loc=rfig, 
                  caption.font.size=my.caption.font.size,
                  foot.font.size=my.foot.font.size,
                  title="Variable importance for Major Adverse Cardiovascular Events (MACE) by 5 years of follow-up", 
                  bkm="vimp_mace", 
                  width=8*0.75, height=4.5*0.75)

update_rtf_fig_numbering(paste0(rfig, 'vimp_mace.rtf'), overwrite=T)

################################################################################
# ASSESS VARIABLE IMPORTANCE: Diabetes
################################################################################

# Extract variable names used to fit crf
X_diab <- ite_diab %>% select(-all_of(c('patid', 'ite_diab', 'se_diab')))

# Calculate variable importance
vimp_diab <- data.frame(covar=colnames(X_diab),
                        importance=variable_importance(crf_diab)) 

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

# Need to clean up variable labels and values 
# This splits strings, applies regex to cleanup, and formats numeric values
vimp_diab_label <- vimp_diab %>%
  # Variable = first str before .
  # Values   = second, third, fourth str split by .
  mutate(var =str_split_i(covar, '[.]', 1),
         val1=str_split_i(covar, '[.]', 2),
         val2=str_split_i(covar, '[.]', 3),
         val3=str_split_i(covar, '[.]', 4)) %>%
  # Replace empty splits ("NA") with empty str
  mutate_all(~replace_na(., '')) %>%
  # Combine val1, val2, val3 into val and apply regex to get nice value labels
  mutate(val=paste(val1, val2, val3),
         # replace GTE and LT with appropriate symbols
         val=gsub("GTE", "\u2265", val),
         val=gsub("LT", "<", val),
         # add space before capital letter to separate compounded words
         val=trimws(gsub("([A-Z])", " \\1", val)),
         # add hyphen between two numbers
         val=gsub("(\\d) (\\d)", "\\1-\\2", val),
         # add + to standalone number
         val=ifelse(!grepl('-|\u2265|<', val), gsub("(\\d)$", "\\1+", val), val),
         # replace "and" and duplicate spaces with 1 space
         val=gsub("and", " ", val),
         val=gsub("  ", " ", val),
         # replace 'Y' with 'Yes'
         val=ifelse(val=='Y', 'Yes', val),
         val=ifelse(val=='N', 'No', val)) %>%
  left_join(label_df, by=join_by(var==short))  %>%
  # Remove ugly vals and no longer needed covar column
  select(label, var, val, importance) 

vimp_diab_report <- vimp_diab_label %>%
  group_by(label) %>%
  reframe(var=var, val='Total', importance=mean(importance)) %>%
  bind_rows(vimp_diab_label) %>%
  arrange(label) %>%
  distinct()

vimp_diab_plot <- vimp_diab_report %>%
  filter(val=='Total') %>%
  ggplot(aes(y=reorder(label, importance), x=importance)) +
  geom_col(alpha=0.75, aes(fill=label)) +
  geom_line(group=1) +
  geom_point(group=1) +
  labs(x='Variable importance', y=NULL) +
  guides(fill="none") +
  theme_bw() +
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=16, face='bold'),
        axis.title.y=element_text(size=16, face='bold'),
        axis.text.y=element_text(size=16),
        legend.title=element_text(size=12, face='bold'),
        legend.text=element_text(size=12),
        legend.box.background=element_rect(color='black', linewidth=1),
        plot.margin=unit(c(1, 0.5, 5.5, 5.5), 'pt'))

vimp_diab_plot  

# Save to output directory
ggsave(paste0(rfig, 'f_var_imp/vimp_diab.png'), plot=vimp_diab_plot, width=8, height=4.5)

insertPicture_rtf(plotin=paste0(rfig, 'f_var_imp/vimp_diab.png'), 
                  path=this.path(),
                  fileout_loc=rfig, 
                  caption.font.size=my.caption.font.size,
                  foot.font.size=my.foot.font.size,
                  title="Variable importance for Diabetes Mellitus (DM) by 5 years of follow-up", 
                  bkm="vimp_diab", 
                  width=8*0.75, height=4.5*0.75)


update_rtf_fig_numbering(paste0(rfig, 'vimp_diab.rtf'), overwrite=T)