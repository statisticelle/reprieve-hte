################################################################################
# PROD PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/production/f_kmcurves.r
# WORK PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/working/f_kmcurves.r
# 
# RISK ASSESSMENT: Low
# DOCUMENT ID:     PROG.10084.v3
# 
# PURPOSE:         Create Kaplan-Meier curves for events and censoring
# 
# SOURCE PRGM:     none
# INPUT:           
# OUTPUT:          /home/reprieve/POOLED/final/MACE_hte/output/f_kmcurves/km_mace.png
#                                                                         km_diab.png
#                                                                         km_cens.png
#                                                                         km_mace_croi.png
#                  /home/reprieve/POOLED/final/MACE_hte/output/km_mace.rtf
#                                                              km_diab.rtf
#                                                              km_cens.rtf
# PACKAGES USED:   tidyverse, survival, officer, table1, this.path, gridExtra
# FUNCTIONS USED:  read.sas, functions/r_word_Helpers.R, functions/rtf_numbering.R
# EXEMPTIONS:      none
# 
# AUTHOR:          Emma Davies Smith
# CREATION DATE:   21NOV2025
# 
# NOTES:           ** Need to figure out how to insert pngs into rtfs for report
# MODIFICATIONS:   none 
################################################################################

################################################################################
# SETUP
################################################################################

# Project setupfile
source("/home/reprieve/POOLED/final/MACE_hte/programs/setupfile.R")

# Libraries
library(tidyverse)                # data manipulation and plotting
library(survival)                 # survival analysis
library(officer   ,lib=rwgdir)
library(table1    ,lib=rwgdir)
library(this.path ,lib=rwgdir)
library(gridExtra, lib=pkgdir)

# Reporting helper functions (from RWG)
source("/home/reprieve/POOLED/final/MACE_hte/functions/r_word_helpers.R")
source("/home/reprieve/POOLED/final/MACE_hte/functions/rtf_numbering.R")

# CBAR's read.sas function
source("/home/programs/R/read.sas.R")

################################################################################
# LOAD DATA 
################################################################################

# Load dataset: adgrf
adgrf <- readRDS(paste0(rdata, 'adgrf.RDS'))

# Subset to complete cases
adgrfcc <- adgrf %>% filter(ccfl)

# Time to MACE
adttmace <- read.sas(libname=sdata, dsname="m061_adttmace")

# MACE Components
adce <- read.sas(libname=sdata, dsname="m06_adce")

################################################################################
# FIT + PLOT SURVIVAL CURVES FOR MACE
################################################################################

# Fit KM curves for MACE by randomized arm, among all participants
fitm <- survfit(Surv(avalm, evntm) ~ rxl, data=adgrfcc)
sfitm <- summary(fitm)

fitm_df <- data.frame(x=sfitm$time, 
                      y=sfitm$surv, 
                      l=sfitm$lower,
                      u=sfitm$upper,
                      w=sfitm$strata,
                      r='All participants') 

# Create cumulative incidence (1-KM) plot with confidence intervals
km_mace <- fitm_df %>% 
  mutate(w = ifelse(w=='rxl=Pitavastatin', 'Pitavastatin', 'Placebo')) %>%
  ggplot(aes(x=x, y=1-y, col=w)) + 
  geom_step(linewidth=0.75) + 
  geom_ribbon(aes(ymin=1-l, ymax=1-u, fill=w), alpha=0.25, col=NA) +
  theme_bw() +
  theme(legend.position='bottom', legend.direction='horizontal') +
  labs(y='Cumulative MACE incidence', x='Follow-up year', 
       col='Arm', fill='Arm') +
  theme(axis.text.x=element_text(size=10), 
        strip.text=element_text(size=12, face='bold'),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=12, face='bold'),
        axis.title.y=element_text(size=12, face='bold'),
        legend.title=element_text(size=10, face='bold'),
        legend.text=element_text(size=10),
        legend.box.background=element_rect(color='black', linewidth=1)) +
  coord_cartesian(xlim=c(0,5), ylim=c(0, 0.10)) +
  scale_y_continuous(labels=function(x) paste0(x*100, '%'),
                     limits=c(0, 0.11), breaks=seq(0, 0.10, 0.025)) +
  scale_x_continuous(limits=c(0, 7), breaks=seq(0, 6, 1)) 

km_mace

# Save to output directory
ggsave(paste0(rfig, 'f_kmcurves/km_mace.png'), plot=km_mace, width=8, height=5)

insertPicture_rtf(plotin=paste0(rfig, 'f_kmcurves/km_mace.png'), 
                   path=this.path(),
                   fileout_loc=rfig, 
                   caption.font.size=my.caption.font.size,
                   foot.font.size=my.foot.font.size,
                   title="Cumulative incidence of MACE", 
                   footnote="Shading represents 95% confidence interval.",
                   bkm="km_mace", 
                   width=8*0.75, height=5*0.75)

update_rtf_fig_numbering(paste0(rfig, 'km_mace.rtf'), overwrite=T)

################################################################################
# SURVIVAL PROBABILITIES: MACE
################################################################################

print("Yearly survival probabilities for MACE among all participants")
summary(fitm, times=seq(0,6,1))

################################################################################
# FIT + PLOT SURVIVAL CURVES FOR DIABETES
################################################################################

fitd <- survfit(Surv(avald, evntd) ~ rxl, data=adgrfcc)
sfitd <- summary(fitd)

fitd_df <- data.frame(x=sfitd$time, 
                      y=sfitd$surv, 
                      l=sfitd$lower,
                      u=sfitd$upper,
                      w=sfitd$strata,
                      r='All participants') 

km_diab <- fitd_df %>% 
  mutate(w = ifelse(w=='rxl=Pitavastatin', 'Pitavastatin', 'Placebo')) %>%
  ggplot(aes(x=x, y=1-y, col=w)) + 
  geom_step(linewidth=0.75) + 
  geom_ribbon(aes(ymin=1-l, ymax=1-u, fill=w), alpha=0.25, col=NA) +
  theme_bw() +
  theme(legend.position='bottom', legend.direction='horizontal') +
  labs(y='Cumulative DM incidence', x='Follow-up year', 
       col='Arm', fill='Arm') +
  theme(axis.text.x=element_text(size=10), 
        strip.text=element_text(size=12, face='bold'),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=12, face='bold'),
        axis.title.y=element_text(size=12, face='bold'),
        legend.title=element_text(size=10, face='bold'),
        legend.text=element_text(size=10),
        legend.box.background=element_rect(color='black', linewidth=1)) +
  coord_cartesian(xlim=c(0,5), ylim=c(0, 0.10)) +
  scale_y_continuous(labels=function(x) paste0(x*100, '%'),
                     limits=c(0, 0.10), breaks=seq(0, 0.10, 0.025)) +
  scale_x_continuous(limits=c(0, 7), breaks=seq(0, 6, 1)) 

km_diab

ggsave(paste0(rfig, 'f_kmcurves/km_diab.png'), plot=km_diab, width=8, height=5)

insertPicture_rtf(plotin=paste0(rfig, 'f_kmcurves/km_diab.png'), 
                   path=this.path(),
                   fileout_loc=rfig, 
                   caption.font.size=my.caption.font.size,
                   foot.font.size=my.foot.font.size,
                   title="Cumulative incidence of Diabetes Mellitus (DM)",
                   footnote="Shading represents 95% confidence interval.",
                   bkm="km_diab", 
                   width=8*0.75, height=5*0.75)

update_rtf_fig_numbering(paste0(rfig, 'km_diab.rtf'), overwrite=T)

################################################################################
# SURVIVAL PROBABILITIES: DIABETES
################################################################################

print("Yearly survival probabilities for diabetes")
summary(fitd, times=seq(0,6,1))

################################################################################
# COMBINED PLOT
################################################################################

km_grid <- grid.arrange(km_mace, km_diab, ncol=2)
ggsave(paste0(rfig, 'f_kmcurves/km_grid.png'), plot=km_grid, width=12, height=6)
insertPicture_rtf(plotin=paste0(rfig, 'f_kmcurves/km_grid.png'), 
                  path=this.path(),
                  fileout_loc=rfig, 
                  caption.font.size=my.caption.font.size,
                  foot.font.size=my.foot.font.size,
                  title="Cumulative incidence of MACE and Diabetes Mellitus (DM)",
                  footnote="Shading represents 95% confidence interval.",
                  bkm="km", 
                  width=12*0.5, height=6*0.5)

update_rtf_fig_numbering(paste0(rfig, 'km.rtf'), overwrite=T)

################################################################################
# FIT + PLOT SURVIVAL CURVES FOR CENSORING
################################################################################

fitc <- survfit(Surv(avalm, 1-evntm) ~ rxl, data=adgrfcc)
sfitc <- summary(fitc)

fitc_df <- data.frame(x=sfitc$time, 
                      y=sfitc$surv, 
                      l=sfitc$lower,
                      u=sfitc$upper,
                      w=sfitc$strata,
                      r='All participants') 

km_cens <- fitc_df %>% 
  mutate(w = ifelse(w=='rxl=Pitavastatin', 'Pitavastatin', 'Placebo')) %>%
  ggplot(aes(x=x, y=1-y, col=w)) + 
  geom_step(linewidth=0.75) + 
  geom_ribbon(aes(ymin=1-l, ymax=1-u, fill=w), alpha=0.25, col=NA) +
  theme_bw() +
  theme(legend.position=c(0.15, 0.75)) +
  labs(y='Cumulative incidence of censoring', x='Follow-up year', 
       col='Randomized arm', fill='Randomized arm') +
  theme(axis.text.x=  element_text(size=12), 
        strip.text=   element_text(size=12, face='bold'),
        axis.text.y=  element_text(size=12),
        axis.title.x= element_text(size=12, face='bold'),
        axis.title.y= element_text(size=12, face='bold'),
        legend.title= element_text(size=12, face='bold'),
        legend.text=  element_text(size=12),
        legend.box.background=element_rect(color='black', linewidth=1),
        plot.margin=unit(c(5.5, 0.5, 5.5, 0.5), 'pt')) +
  coord_cartesian(xlim=c(0,6), ylim=c(0, 0.80)) +
  scale_y_continuous(labels=function(x) paste0(x*100, '%'),
                     limits=c(0, 0.80), breaks=seq(0, 0.80, 0.10)) +
  scale_x_continuous(limits=c(0, 7), breaks=seq(0, 6, 1)) 

km_cens

ggsave(paste0(rfig, 'f_kmcurves/km_cens.png'), plot=km_cens, width=8, height=3.25)

insertPicture_rtf(plotin=paste0(rfig, 'f_kmcurves/km_cens.png'), 
                   path=this.path(),
                   fileout_loc=rfig, 
                   caption.font.size=my.caption.font.size,
                   foot.font.size=my.foot.font.size,
                   title="Cumulative incidence of censoring. Shading represents 95% confidence interval.", 
                   bkm="km_cens", 
                   width=8*0.75, height=3.25*0.75)

update_rtf_fig_numbering(paste0(rfig, 'km_cens.rtf'), overwrite=T)

################################################################################
# SURVIVAL PROBABILITIES: DIABETES
################################################################################

print("Yearly survival probabilities for censoring")
summary(fitc, times=seq(0,6,1))

################################################################################
# EVENT COUNTS
################################################################################

adgrfcc <- adgrfcc %>%
  mutate(evntm5=case_when(evntm==1 & avalm <= 5 ~ 1,
                          avalm >= 5 ~ 0,
                          evntm==0 & avalm < 5 ~ NA)) 

print('MACE events within 5 years, by arm')
adgrfcc %>%
  group_by(rxl, evntm5) %>%
  summarize(n())

adgrfcc <- adgrfcc %>%
  mutate(evntd5=case_when(evntd==1 & avald <= 5 ~ 1,
                          avald >= 5 ~ 0,
                          evntd==0 & avald < 5 ~ NA))

print('Diabetes events within 5 years, by arm')
adgrfcc %>%
  group_by(rxl, evntd5) %>%
  summarize(n())

print('MACE and Diabetes events within 5 years, by arm') 
adgrfcc %>%
  group_by(rxl, evntd5, evntm5) %>%
  summarize(n())

print('Diabetes before MACE within 5 years, by arm')
adgrfcc %>%
  filter(evntd5==1 & evntm5==1 & avald < avalm) %>%
  group_by(rxl) %>%
  summarize(n())
  

################################################################################
# MACE COMPONENTS
################################################################################

grfce <- adce %>%
  filter(patid %in% adgrfcc$patid) %>%
  filter(AOCCPMFL=='Y') %>%
  left_join(adgrfcc)

print('MACE component events within 5 years, by arm')
grfce %>%
  filter(evntm5 == 1) %>%
  group_by(rxl, cecat, class) %>%
  summarize(n=n())%>%
  arrange(rxl, -n) %>%
  print(n=30)

print('MACE component events within 5 years, overall')
grfce %>%
  filter(evntm5 == 1) %>%
  mutate(class=ifelse(cecat=='Death', 'Death', class)) %>%
  group_by(cecat, class) %>%
  summarize(n=n())%>%
  ungroup() %>%
  mutate(p=n/sum(n)) %>%
  arrange(-n) 
