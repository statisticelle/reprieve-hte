################################################################################
# PROD PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/f_quadrants.r
# WORK PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/working/f_quadrants.r
# 
# RISK ASSESSMENT: LOW
# DOCUMENT ID:     PROG.10084.v3
# 
# PURPOSE:         Visualize risk-benefit quadrants, according to combinations of key risk factors
# 
# SOURCE PRGM:     none
# INPUT:           /sas_datasets/Rderived/adgrf.RDS, ite_mace.RDS, ite_diab.RDS
# OUTPUT:          /output/figures/f_quadrants/quadrants.png
#                  /output/figures/quadrants.rtf
#
# PACKAGES USED:   none
# FUNCTIONS USED:  none
# EXEMPTIONS:      none
# 
# AUTHOR:          Emma Davies Smith
# CREATION DATE:   10DEC2025
# 
# NOTES:           none
# MODIFICATIONS:   none 
################################################################################

################################################################################
# SETUP
################################################################################

# Project setupfile
source("/home/reprieve/POOLED/final/MACE_hte/programs/setupfile.R")

# Libraries
library(tidyverse)                # data manipulation and plotting
library(officer   ,lib=rwgdir)
library(table1    ,lib=rwgdir)
library(this.path ,lib=rwgdir)
library(flextable ,lib=rwgdir)

# Reporting helper functions (from RWG)
source("/home/reprieve/POOLED/final/MACE_hte/functions/r_word_helpers.R")
source("/home/reprieve/POOLED/final/MACE_hte/functions/rtf_numbering.R")

################################################################################
# LOAD DATA 
################################################################################

# Load analysis dataset
adgrf <- readRDS(paste0(rdata, 'adgrf.RDS'))

# Load ITE predictions
ite_mace <- readRDS(paste0(rdata, 'ite_mace.RDS')) 
ite_diab <- readRDS(paste0(rdata, 'ite_diab.RDS'))

# Combine ITE predictions
ite <- ite_mace %>% select(patid, ite_mace, se_mace) %>%
  left_join(ite_diab %>% select(patid, ite_diab, se_diab), by='patid') %>%
  filter(!is.na(ite_diab)) %>%
  # create prediction interval limits for ites (for exploration)
  mutate(lower_mace = ite_mace - 1.96 * se_mace,
         upper_mace = ite_mace + 1.96 * se_mace,
         lower_diab = ite_diab - 1.96 * se_diab,
         upper_diab = ite_diab + 1.96 * se_diab) %>%
  # add ascvd < 5% boolean and covariates
  left_join(adgrf %>% select(patid, egfrfl, fglucfl, bmigp6), by='patid')

################################################################################
# Compute quadrant probabilities (simplifying assumption: ITEs are indepdnent)
################################################################################

pr_quad <- function(quad, ite_diab=ite_diab, se_diab=se_diab, ite_mace=ite_mace, se_mace=se_mace){
  if (quad=='tl'){
    # PR(MACE ITE < 0) AND PR(DIAB ITE > 0)
    pnorm(0, mean=ite_mace, sd=se_mace) * (1-pnorm(0, mean=ite_diab, sd=se_diab))
  } else if (quad=='tr') {
    # PR (MACE ITE > 0) AND PR(DIAB ITE > 0)
    (1-pnorm(0, mean=ite_diab, sd=se_diab)) * (1-pnorm(0, mean=ite_mace, sd=se_mace))
  } else if (quad=='br') {
    # PR (MACE ITE > 0) AND PR(DIAB ITE < 0)
    (1-pnorm(0, mean=ite_mace, sd=se_mace)) * (pnorm(0, mean=ite_diab, sd=se_diab))
  } else if (quad=='bl') {
    # PR (MACE ITE < 0) AND PR(DIAB ITE < 0)
    (pnorm(0, mean=ite_diab, sd=se_diab)) * (pnorm(0, mean=ite_mace, sd=se_mace))
  } else {
    print('Specify valid quadrant')
  }
}

quaddf <- ite %>%
  mutate(obs_quad = case_when(ite_mace < 0 & ite_diab > 0 ~ 'tl',
                              ite_mace > 0 & ite_diab > 0 ~ 'tr',
                              ite_mace > 0 & ite_diab < 0 ~ 'br',
                              ite_mace < 0 & ite_diab < 0 ~ 'bl')) %>%
  mutate(pr_tl = pr_quad(quad='tl', ite_diab, se_diab, ite_mace, se_mace),
         pr_tr = pr_quad(quad='tr', ite_diab, se_diab, ite_mace, se_mace),
         pr_br = pr_quad(quad='br', ite_diab, se_diab, ite_mace, se_mace),
         pr_bl = pr_quad(quad='bl', ite_diab, se_diab, ite_mace, se_mace),
         pr_quad = case_when(obs_quad == 'tl' ~ pr_tl,
                             obs_quad == 'tr' ~ pr_tr,
                             obs_quad == 'br' ~ pr_br,
                             obs_quad == 'bl' ~ pr_bl)) %>%
  mutate(pr_tier = cut(pr_quad, breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c('<25%', '25-50%', '50-75%', '>75%'))) %>%
  mutate(egfrfl =ifelse(egfrfl=='LT90', 1, 0),
         gbfl   =ifelse(fglucfl=='GTE100' | bmigp6=='Obese', 1, 0)) %>%
  mutate(label  = case_when(egfrfl==1 & gbfl==0 ~ 'Lower eGFR (25%)',
                            egfrfl==0 & gbfl==1 ~ 'Higher glucose or obese (19%)',
                            egfrfl==1 & gbfl==1 ~ 'Both (16%)',
                            egfrfl==0 & gbfl==0 ~ 'Neither (41%)')) %>%
  mutate(label = factor(label, levels=c('Neither (41%)', 'Lower eGFR (25%)', 'Higher glucose or obese (19%)', 'Both (16%)')))

################################################################################
# Tabulate quadrant frequencies
################################################################################

print('Quadrant frequencies:')
table(quaddf$obs_quad)
round(table(quaddf$obs_quad)/sum(table(quaddf$obs_quad)),2)*100

print('Cluster frequencies:')
table(quaddf$label)
round(table(quaddf$label)/sum(table(quaddf$label)),2)*100

quaddf %>%
  group_by(label) %>%
  summarize(sprintf('%.1f', mean(ite_mace)*100),
            sprintf('%.1f', mean(ite_diab)*100))
  

################################################################################
# Risk-benefit matrix with risk factors and with bivariate CIs
################################################################################

colors <- c('Neither (41%)'='#7CAE00', 'Lower eGFR (25%)'='#00BFC4', 'Higher glucose or obese (19%)'='#F8766D', 'Both (16%)'='#C77CFF')

quadggci <- ggplot() +
  geom_point(size=3, aes(y=ite_diab, x=ite_mace, col='Neither (41%)', alpha=0.75), data=quaddf[quaddf$label=='Neither (41%)',]) +
  geom_point(size=3, aes(y=ite_diab, x=ite_mace, col='Lower eGFR (25%)', alpha=0.75), data=quaddf[quaddf$label=='Lower eGFR (25%)',]) +
  geom_point(size=3, aes(y=ite_diab, x=ite_mace, col='Higher glucose or obese (19%)', alpha=0.75), data=quaddf[quaddf$label=='Higher glucose or obese (19%)',]) +
  geom_point(size=3, aes(y=ite_diab, x=ite_mace, col='Both (16%)', alpha=0.75), data=quaddf[quaddf$label=='Both (16%)',]) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  stat_ellipse(aes(y=ite_diab, x=ite_mace, col=label, fill=label), type='t', geom='polygon', level=0.95, alpha=0.15, linewidth=1, data=quaddf) +
  labs(y='Estimated CATE for diabetes', x='Estimated CATE for MACE', col='Key risk factors') +
  scale_y_continuous(labels=function(x) paste0(x*100, '%'),
                     limits=c(-0.005, 0.0425), breaks=seq(-0.01, 0.05, 0.01)) +
  scale_x_continuous(labels=function(x) paste0(x*100, '%'),
                     limits=c(-0.0425, 0.005), breaks=seq(-0.05, 0.01, 0.01)) +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
  guides(fill='none', alpha='none') +
  theme_bw() +
  theme(legend.position='inside',
        legend.position.inside=c(0.25, 0.85),
        legend.box='horizontal',
        legend.background=element_rect(color='black', linetype='solid'),
        axis.text.x=element_text(size=16), 
        strip.text=element_text(size=16, face='bold'),
        axis.text.y=element_text(size=16),
        axis.title.x=element_text(size=16, face='bold'),
        axis.title.y=element_text(size=16, face='bold'),
        legend.title=element_text(size=16, face='bold'),
        legend.text=element_text(size=16),
        legend.spacing = unit(0.1, 'cm'),
        legend.margin=margin(0.1,0.1,0.1,0.1, unit='cm'),
        legend.key.spacing=unit(0.05, 'cm')) +
  annotate('label', x=-0.035, y= 0, fontface='bold',  size=16/.pt, label='Decreased MACE risk\nwith pitavastatin') +
  annotate('label', x=-0.005, y= 0.0425, fontface='bold',  size=16/.pt, label='Increased diabetes risk') 

quadggci
# Save ggplot as png and output
ggsave('/home/reprieve/POOLED/final/MACE_hte/publication/AIDS/output/fig2.png', plot=quadggci, width=10, height=8)


