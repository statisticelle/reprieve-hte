################################################################################
# PROD PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/production/f_ite_dist.r
# WORK PROGRAM:    /home/reprieve/POOLED/final/MACE_hte/programs/tables_figures/working/f_dist.r
# 
# RISK ASSESSMENT: Low
# DOCUMENT ID:     PROG.10084.v3
# 
# PURPOSE:         Visualize distributions of ITEs for both events
# 
# SOURCE PRGM:     none
# INPUT:           rdata/adgrf.RDS
#                        ite_mace.RDS
#                        ite_diab.RDS
# OUTPUT:          rfig/f_ite_boxplot/ite_density.png
#                                     box_mace.png
#                                     box_diab.png
#                  rfig/ite_density.rtf
#                       box_mace.rtf
#                       box_diab.rtf
#
# PACKAGES USED:   none
# FUNCTIONS USED:  none
# EXEMPTIONS:      none
# 
# AUTHOR:          Emma Davies Smith
# CREATION DATE:   02DEC2025
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
library(gridExtra, lib=pkgdir)
library(officer   ,lib=rwgdir)
library(table1    ,lib=rwgdir)
library(this.path ,lib=rwgdir)

# Reporting helper functions (from RWG)
source("/home/reprieve/POOLED/final/MACE_hte/functions/r_word_helpers.R")
source("/home/reprieve/POOLED/final/MACE_hte/functions/rtf_numbering.R")

################################################################################
# LOAD DATA 
################################################################################

# Load analysis dataset
adgrf <- readRDS(paste0(rdata, 'adgrf.RDS'))
adgrfcc <- adgrf %>% filter(ccfl)


# Load ITE predictions
ite_mace <- readRDS(paste0(rdata, 'ite_mace.RDS')) 
ite_diab <- readRDS(paste0(rdata, 'ite_diab.RDS'))

# Combine ITE predictions
ite <- ite_mace %>% select(patid, ite_mace, se_mace) %>%
  left_join(ite_diab %>% select(patid, ite_diab, se_diab), by='patid') %>%
  # add ascvd < 5% boolean, treatment arm, and important covariates
  left_join(adgrf %>% select(patid, rxl, gbdrace, bmigp6, egfrfl), by='patid')

# Color palette
pal <- c(
  "All participants"  = "gray",
  "Both"              = "#F8766D",
  "Higher Glucose"    = "#E68613",
  "Obese"             = "#E68613",
  "Neither"           = "#7CAE00",
  "Higher eGFR"       = "#7CAE00",
  "Lower eGFR"        = "#F8766D"
)

################################################################################
# BOXPLOTS OF ITES: MACE
###############################################################################

# Among all participants -------------------------------------------------------

all_mace_df <- ite %>%
  mutate(class='All participants') 

lab_all_mace <- all_mace_df %>%
  group_by(class) %>%
  summarize(n=n(),
            min=quantile(ite_mace, p=0.025)*100,
            q1=quantile(ite_mace, p=0.25)*100,
            q2=quantile(ite_mace, p=0.50)*100,
            q3=quantile(ite_mace, p=0.75)*100,
            max=quantile(ite_mace, p=0.975)*100) %>%
  mutate(pct=ifelse(n/sum(n) < 0.01, round(n/sum(n)*100, 1), round(n/sum(n)*100, 0)),
         ntxt=paste0(n, ' (', pct, '%)'),
         iqrtxt=paste0(sprintf('%.1f', q1), '% to ', sprintf('%.1f', q3), '%'),
         q2txt=paste0(sprintf('%.1f', q2), '%'),
         rngtxt=paste0(sprintf('%.1f', min), '% to ', sprintf('%.1f', max), '%'),
         txt=paste0(class, '\n', ntxt, '\n', q2txt, '\n', iqrtxt, '\n', rngtxt)) %>%
  select(class, txt)

lab_all_mace <- rbind(lab_all_mace,
                      c(class='', txt='\nN (%)\nMedian\nQ1 to Q3\n5th to 95th %tile'))

box_mace_all <- all_mace_df %>%
  ggplot(aes(y=ite_mace, x=class, fill=class)) +
  geom_hline(yintercept=0, alpha=0.5, linetype=2) +
  geom_boxplot(outlier.alpha=0.25) +
  geom_label(data=lab_all_mace, size=14/.pt, fill='white', alpha=0.85, label.size=0,
             aes(label=txt, y=0.025, x=class, group=class, 
                 fontface=ifelse(class=='', 'bold', 'plain'))) +
  geom_segment(aes(x='', xend='', y=0, yend=-0.05), arrow=arrow(length=unit(0.3, "in"))) +
  annotate("text", x="", y=-0.02, hjust=0.5, vjust=-1, label="Decreasing", angle=90, size=14/.pt) +
  annotate("text", x="", y=-0.02, hjust=0.5, vjust=1.5,  label="MACE risk", angle=90, size=14/.pt, fontface='italic') +
  scale_x_discrete(drop=F, labels=NULL) +
  scale_y_continuous(labels=function(x) paste0(x*100, '%'),
                     limits=c(-0.05, 0.05), breaks=c(-0.05, 0, 0.05)) +
  scale_fill_manual(values=pal) +
  labs(x=NULL, y='CATE for 5-year MACE', fill=NULL) +
  guides(fill='none') +
  theme_bw() +
  theme(axis.text.x=element_text(size=16, face='bold'),
        axis.text.y=element_text(size=16, colour='black'),
        axis.title.x=element_text(size=16, face='bold'),
        axis.title.y=element_text(size=16, face='bold'),
        legend.title=element_text(size=12, face='bold'),
        legend.text=element_text(size=12),
        legend.box.background=element_rect(color='black', linewidth=1),
        plot.margin=unit(c(1, 0.5, 5.5, 5.5), 'pt')) 

box_mace_all 

# By risk factor ---------------------------------------------------------------

box_mace_df <- ite %>%
  left_join(adgrfcc %>% select(patid, egfrfl)) %>%
  mutate(class=case_when(egfrfl=='LT90' ~ 'Lower eGFR',
                         egfrfl=='GTE90' ~ 'Higher eGFR')) 

lab_mace_df <- box_mace_df %>%
  group_by(class) %>%
  summarize(n=n(),
            min=quantile(ite_mace, p=0.025)*100,
            q1=quantile(ite_mace, p=0.25)*100,
            q2=quantile(ite_mace, p=0.50)*100,
            q3=quantile(ite_mace, p=0.75)*100,
            max=quantile(ite_mace, p=0.975)*100) %>%
  mutate(pct=ifelse(n/sum(n) < 0.01, round(n/sum(n)*100, 1), round(n/sum(n)*100, 0)),
         ntxt=paste0(n, ' (', pct, '%)'),
         iqrtxt=paste0(sprintf('%.1f', q1), '% to ', sprintf('%.1f', q3), '%'),
         q2txt=paste0(sprintf('%.1f', q2), '%'),
         rngtxt=paste0(sprintf('%.1f', min), '% to ', sprintf('%.1f', max), '%'),
         txt=paste0(class, '\n', ntxt, '\n', q2txt, '\n', iqrtxt, '\n', rngtxt)) %>%
  select(class, txt)

box_mace <- box_mace_df %>%
  ggplot(aes(y=ite_mace, x=class, fill=class)) +
  geom_hline(yintercept=0, alpha=0.5, linetype=2) +
  geom_boxplot(outlier.alpha=0.25) +
  geom_label(data=lab_mace_df, size=14/.pt, fill='white', alpha=0.85, label.size=0,
             aes(label=txt, y=0.025, x=class, group=class, 
                 fontface=ifelse(class=='', 'bold', 'plain'))) +
  scale_x_discrete(drop=F, labels=NULL) +
  scale_y_continuous(labels=NULL,
                     limits=c(-0.05, 0.05), breaks=c(-0.05, 0, 0.05)) +
  scale_fill_manual(values=pal) +
  labs(x=NULL, y=NULL, fill=NULL) +
  guides(fill='none') +
  theme_bw() +
  theme(axis.text.x=element_text(size=16, face='bold'),
        axis.text.y=element_text(size=16, colour='black'),
        axis.title.x=element_text(size=16, face='bold'),
        axis.title.y=element_text(size=16, face='bold'),
        legend.title=element_text(size=12, face='bold'),
        legend.text=element_text(size=12),
        legend.box.background=element_rect(color='black', linewidth=1),
        plot.margin=unit(c(1, 0.5, 5.5, 5.5), 'pt')) 

box_mace

# Combine  -----
mace_grid <- grid.arrange(box_mace_all, box_mace, ncol=2, widths=c(1,2))


################################################################################
# BOXPLOTS OF ITES: Diabetes
################################################################################

summary(ite$ite_diab)*100
table(ite$ite_diab > 0) 

# Among all participants -------------------------------------------------------

all_diab_df <- ite %>%
  mutate(class='All participants') 

lab_all_diab <- all_diab_df %>%
  group_by(class) %>%
  summarize(n=n(),
            min=quantile(ite_diab, p=0.025)*100,
            q1=quantile(ite_diab, p=0.25)*100,
            q2=quantile(ite_diab, p=0.50)*100,
            q3=quantile(ite_diab, p=0.75)*100,
            max=quantile(ite_diab, p=0.975)*100) %>%
  mutate(pct=ifelse(n/sum(n) < 0.01, round(n/sum(n)*100, 1), round(n/sum(n)*100, 0)),
         ntxt=paste0(n, ' (', pct, '%)'),
         iqrtxt=paste0(sprintf('%.1f', q1), '% to ', sprintf('%.1f', q3), '%'),
         q2txt=paste0(sprintf('%.1f', q2), '%'),
         rngtxt=paste0(sprintf('%.1f', min), '% to ', sprintf('%.1f', max), '%'),
         txt=paste0(class, '\n', ntxt, '\n', q2txt, '\n', iqrtxt, '\n', rngtxt)) %>%
  select(class, txt)


lab_all_diab <- rbind(lab_all_diab,
                      c(class='', txt='\nN (%)\nMedian\nQ1 to Q3\n5th to 95th %tile'))

box_diab_all <- all_diab_df %>%
  ggplot(aes(y=ite_diab, x=class, fill=class)) +
  geom_hline(yintercept=0, alpha=0.5, linetype=2) +
  geom_boxplot(outlier.alpha=0.25) +
  geom_label(data=lab_all_diab, size=14/.pt, fill='white', alpha=0.85, label.size=0,
             aes(label=txt, y=-0.025, x=class, group=class, 
                 fontface=ifelse(class=='', 'bold', 'plain'))) +
  geom_segment(aes(x='', xend='', y=0, yend=0.05), arrow=arrow(length=unit(0.3, "in"))) +
  annotate("text", x="", y=0.025, hjust=0.5, vjust=-1, label="Increasing", angle=90, size=14/.pt) +
  annotate("text", x="", y=0.025, hjust=0.5, vjust=1.5, label="diabetes risk", angle=90, size=14/.pt, fontface='italic') +
  scale_x_discrete(drop=F, labels=NULL) +
  scale_y_continuous(labels=function(x) paste0(x*100, '%'),
                     limits=c(-0.05, 0.05), breaks=c(-0.15, -0.10, -0.05, 0, 0.05, 0.10)) +
  scale_fill_manual(values=pal) +
  labs(x=NULL, y='CATE for 5-year diabetes', fill=NULL) +
  guides(fill='none') +
  theme_bw() +
  theme(axis.text.x=element_text(size=16, face='bold'),
        axis.text.y=element_text(size=16, colour='black'),
        axis.title.x=element_text(size=16, face='bold'),
        axis.title.y=element_text(size=16, face='bold'),
        legend.title=element_text(size=12, face='bold'),
        legend.text=element_text(size=12),
        legend.box.background=element_rect(color='black', linewidth=1),
        plot.margin=unit(c(1, 0.5, 5.5, 5.5), 'pt')) 

box_diab_all

# By risk factor ---------------------------------------------------------------
box_diab_df <- ite %>%
  filter(!is.na(ite_diab)) %>%
  left_join(adgrfcc %>% select(patid, fglucfl, bmigp6)) %>%
  mutate(class=case_when(fglucfl=='GTE100' & bmigp6!='Obese' ~ 'Higher Glucose',
                         fglucfl=='GTE100' & bmigp6=='Obese' ~ 'Both',
                         fglucfl=='LT100'  & bmigp6!='Obese' ~ 'Neither',
                         fglucfl=='LT100'  & bmigp6=='Obese' ~ 'Obese',
                         TRUE ~ 'Other')) %>%
  mutate(class=factor(class, levels=c('Both', 'Higher Glucose', 'Obese', 'Neither')))

lab_diab_df <- box_diab_df %>%
  group_by(class) %>%
  summarize(n=n(),
            min=quantile(ite_diab, p=0.025)*100,
            q1=quantile(ite_diab, p=0.25)*100,
            q2=quantile(ite_diab, p=0.50)*100,
            q3=quantile(ite_diab, p=0.75)*100,
            max=quantile(ite_diab, p=0.975)*100) %>%
  mutate(pct=ifelse(n/sum(n) < 0.01, round(n/sum(n)*100, 1), round(n/sum(n)*100, 0)),
         ntxt=paste0(n, ' (', pct, '%)'),
         iqrtxt=paste0(sprintf('%.1f', q1), '% to ', sprintf('%.1f', q3), '%'),
         q2txt=paste0(sprintf('%.1f', q2), '%'),
         rngtxt=paste0(sprintf('%.1f', min), '% to ', sprintf('%.1f', max), '%'),
         txt=paste0(class, '\n', ntxt, '\n', q2txt, '\n', iqrtxt, '\n', rngtxt)) %>%
  select(class, txt)

box_diab <- box_diab_df %>%
  ggplot(aes(y=ite_diab, x=class, fill=class)) +
  geom_hline(yintercept=0, alpha=0.5, linetype=2) +
  geom_boxplot(outlier.alpha=0.25) +
  geom_label(data=lab_diab_df, size=14/.pt, fill='white', alpha=0.85, label.size=0,
             aes(label=txt, y=-0.025, x=class, group=class, 
                 fontface=ifelse(class=='', 'bold', 'plain'))) +
  scale_x_discrete(drop=F, labels=NULL) +
  scale_y_continuous(labels=NULL,
                     limits=c(-0.05, 0.05), breaks=c(-0.15, -0.10, -0.05, 0, 0.05, 0.10)) +
  scale_fill_manual(values=pal) +
  labs(x=NULL, y=NULL, fill=NULL) +
  guides(fill='none') +
  theme_bw() +
  theme(axis.text.x=element_text(size=16, face='bold'),
        axis.text.y=element_text(size=16, colour='black'),
        axis.title.x=element_text(size=16, face='bold'),
        axis.title.y=element_text(size=16, face='bold'),
        legend.title=element_text(size=12, face='bold'),
        legend.text=element_text(size=12),
        legend.box.background=element_rect(color='black', linewidth=1),
        plot.margin=unit(c(1, 0.5, 5.5, 5.5), 'pt')) 

box_diab

diab_grid <- grid.arrange(box_diab_all, box_diab, ncol=2, widths=c(1,2))

# Create figure 1 and save--
fig1 <- grid.arrange(mace_grid, diab_grid, ncol=1)
ggsave(plot=fig1, '/home/reprieve/POOLED/final/MACE_hte/publication/AIDS/output/fig1.png',
       width=14, height=9, dpi="retina")

