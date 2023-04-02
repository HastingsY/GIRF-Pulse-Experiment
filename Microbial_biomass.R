##Microbial Biomass C, N, and C:N
##Yvette Hastings
##April 17, 2022
##Updated 4/2/2023

## Data analysis sections; quick access using shift + alt + J
## 1. Data clean-up
## 2. Microbial Biomass C
## 3. Microbial Biomass N
## 4. Microbial Biomass C:N

##load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(stringr)
library(rstatix)
library(ggpubr)
library(reshape)
library(plyr)
library(datarium)
library(sjPlot)
library(ggpmisc)
library(agricolae)
library(emmeans)


# Data clean-up -----------------------------------------------------------

##load data
all_pulse <- read_excel("Data/2020-2021 Soils Master File.xlsx", sheet = 'Pulse Data Masterfile (R)') 

##calculate C:N
all_pulse$molC <- all_pulse$mb.doc/12.0107
all_pulse$molN <- all_pulse$mb.tdn/14.0067
all_pulse$CN <- all_pulse$molC/all_pulse$molN

##factor sampling dates, sampling campaign, and treatment for graphing
all_pulse$day.of.pulse[all_pulse$day.of.pulse == 'pre-pulse'] <- "Pre-pulse"
all_pulse$day.of.pulse[all_pulse$day.of.pulse == 'Pre-pulse'] <- "0 day (pre-pulse)"

all_pulse$date <- ifelse(all_pulse$day.of.pulse == "0 day (pre-pulse)", 0,
                         ifelse(all_pulse$day.of.pulse == "3 days post-pulse", 3,
                                ifelse(all_pulse$day.of.pulse == "5 days post-pulse", 5,
                                       ifelse(all_pulse$day.of.pulse == "9 days post-pulse", 9,
                                              ifelse(all_pulse$day.of.pulse == "18 days post-pulse", 16,
                                                     ifelse(all_pulse$day.of.pulse == "1 day post-pulse", 1,
                                                            ifelse(all_pulse$day.of.pulse == "2 days post-pulse", 2,
                                                                   ifelse(all_pulse$day.of.pulse == "7 days post-pulse", 7,16))))))))


all_pulse$date <- factor(all_pulse$date, levels = c(0,1,2,3,5,7,9,16))
all_pulse$sampling.campaign <- factor(all_pulse$sampling.campaign, levels = c("GIRF September 2020 Pulse", "June 2021 GIRF Pulse", "June 2021 TM Natural Pulse"))
all_pulse$treatment <- factor(all_pulse$treatment, levels = c("Diverse", "Grass", "Mulch", "Reference"))


##split dataset to individual pulse events
June2021 <- filter(all_pulse, sampling.campaign == "June 2021 GIRF Pulse")  %>%
  drop_na(mb.doc) ##remove rows missing MBC and MBN data
Sept2020 <- filter(all_pulse, sampling.campaign == "GIRF September 2020 Pulse")
TM2021 <- filter(all_pulse, sampling.campaign == "June 2021 TM Natural Pulse")
TM2021_1 <- TM2021[-c(1),] ##remove MBC and CN outlier
all_pulse1 <- all_pulse[-c(46),] #remove outlier in TM DOC for sampling campaign comparison



# Microbial Biomass C -----------------------------------------------------
##Sept 2020
summary <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(mb.doc, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "mb.doc",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "Microbial Biomass C by Treatment over Sept 2020 Pulse Experiment", ylab = 'DOC (mg C/g soil)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(mb.doc)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(mb.doc)
min(normality$p)

ggqqplot(Sept2020, "mb.doc", ggtheme = theme_bw()) +
  facet_grid(.~treatment, labeller = "label_both")
##normality check p>0.05; qqplot shows majority of points are along the regression line; accept normality check 

###ANOVA
SeptDOC.aov <- aov(mb.doc~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SeptDOC.aov)
##No significance from within effects


res.aov <- anova_test(
  data = Sept2020, dv = mb.doc, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##Between effect significance by date and interaction of treatment and date

##post hoc on date effect
Sept2020 %>%
  pairwise_t_test(
    mb.doc ~ date,  
    p.adjust.method = "bonferroni"
  )
##No significanc

##post hoc on interaction terms
Sept2020 %>%
  group_by(date) %>%
  pairwise_t_test(
    mb.doc ~treatment, paired = TRUE,
    p.adjust.method = 'bonferroni')
##No significance 


Sept_DOC <- ggline(Sept2020, x="date", y = "mb.doc", color = "treatment", 
                   add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                   size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  labs(x = '', y= 'Microbial Biomass C (mg C/g soil)') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =150, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=150, label = "Post", color = 'black', size = 9)+ylim (0,150) 



##June 2021
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(mb.doc, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "mb.doc",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "DOC by Treatment over June 2021 Pulse Experiment", ylab = 'DOC (mg C/g soil)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(mb.doc)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(mb.doc)
min(normality$p)

ggqqplot(June2021, "mb.doc", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")
##normality check p>0.05; qqplot shows majority of points are along the regression line; accept normality check 

###ANOVA
JunetDOC.aov <- aov(mb.doc~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JunetDOC.aov)
##No significance within effects

res.aov <- anova_test(
  data = June2021, dv = mb.doc, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance between effects

##June DOC plot
June_DOC <- ggline(June2021, x="date", y = "mb.doc", color = "treatment", 
                   add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                   size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  labs(x = '', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =150, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=150, label = "Post", color = 'black', size = 9) +
  ylim(0,150)


##Todd's Meadow
summary <- TM2021_1 %>%
  group_by(date) %>%
  get_summary_stats(mb.doc, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021_1, x = "date", y = "mb.doc",
          color = "treatment", palette = 'blue', shape = "treatment", add = c("mean_se", "jitter"),
          title = "DOC by Treatment over June TM 2021 Natural Pulse", ylab = 'DOC (mg C/g soil)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021_1 %>%
  group_by(sampling.date) %>%
  identify_outliers(mb.doc)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021_1 %>%
  shapiro_test(mb.doc)
min(normality$p)

ggqqplot(TM2021_1, "mb.doc", ggtheme = theme_bw()) +
  facet_grid(treatment ~ date, labeller = "label_both")
##normality check p>0.05; qqplot shows majority of points are along the regression line; accept normality check 

###ANOVA
TMDOC.aov <- aov(mb.doc~date+ Error(core.id/(date)), data = TM2021_1) ##check interaction terms for sampling date and treatment
summary(TMDOC.aov)
##No significance from within effects

res.aov <- anova_test(
  data = TM2021_1, dv = mb.doc, wid = core.id,
  between = c(date)
)
get_anova_table(res.aov)
##No significance from within effects

##June TM DOC plot
TM_DOC <- ggline(TM2021_1, x="date", y = "mb.doc", color = "treatment", 
                 add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                 size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = 'blue') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  labs(x = '', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =150, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=150, label = "Post", color = 'black', size = 9) + ylim(0,150)


##Between Campaigns
summary <- all_pulse1 %>%
  group_by(sampling.campaign) %>%
  get_summary_stats(mb.doc, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse1, x = "sampling.campaign", y = "mb.doc",
          color = "sampling.campaign", palette = c('grey', 'orange', 'blue'), shape = "sampling.campaign", add = c("mean_se", "jitter"),
          title = "DOC between Experimental and Natural Pulses", ylab = 'DOC (mg C/g soil)', 
          xlab = "", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse1 %>%
  group_by(sampling.campaign) %>%
  identify_outliers(mb.doc)
outliers
##there are two extreme outliers, but were not removed since individual campaign's did not have extreme outliers

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse1 %>%
  group_by(sampling.campaign) %>%
  shapiro_test(mb.doc)
min(normality$p)

ggqqplot(all_pulse1, "mb.doc", ggtheme = theme_bw()) +
  facet_grid(~sampling.campaign, labeller = "label_both")
##majority of points fall along the line

###ANOVA
res.aov <- anova_test(
  data = all_pulse1, dv = mb.doc, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)
##no significance between sampling campaigns 


DOC_sampling.campaign <- ggboxplot(all_pulse1, x = "sampling.campaign", y = "mb.doc",
                                   color = "black", add = c("mean_se", "jitter"),
                                   title = "Microbial Biomass C Between Pulse Events", ylab = 'Microbial Biomass C (mg C/g soil)', 
                                   xlab = "", bxp.errorbar = TRUE,
                                   fill = 'sampling.campaign',
                                   palette =c('grey', 'orange', 'blue'),
                                   legend.title = "Sampling Campaign", legend = "bottom")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 26, face = 'bold', color = 'black'),
        legend.key.size = unit(1, 'cm'), legend.text = element_text(size = 20))+
  ylim(0,200)


# Microbial Biomass N -----------------------------------------------------


##Sept 2020
summary <- Sept2020 %>%
  group_by(sampling.date, treatment) %>%
  get_summary_stats(mb.tdn, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "mb.tdn",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "TDN by Treatment over Sept 2020 Pulse Experiment", ylab = 'TDN (mg N/g soil)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(mb.tdn)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(mb.tdn)
min(normality$p)

ggqqplot(Sept2020, "mb.tdn", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")
##normality check p>0.05; qqplot shows majority of points are along the regression line; accept normality check 

###ANOVA
SeptTDN.aov <- aov(mb.tdn~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SeptTDN.aov)
##No significance from within effects

res.aov <- anova_test(
  data = Sept2020, dv = mb.tdn, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##Between effect significance by date

##post hoc on date effect
Sept2020 %>%
  pairwise_t_test(
    mb.tdn ~ date, 
    p.adjust.method = "bonferroni")

##Sept TDN plot
Sept_TDN <- ggline(Sept2020, x="date", y = "mb.tdn", color = "treatment", title = "",
                   add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                   size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  labs(x = '', y= 'Microbial Biomass N (mg N/g soil)') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =30, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=30, label = "Post", color = 'black', size = 9)+
  geom_bracket(
    xmin = c("0"), xmax = c("16"), 
    y.position = c(15),
    label = c("*"), label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,30)


##June 2021
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(mb.tdn, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "mb.tdn",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "TDN by Treatment over June 2021 Pulse Experiment", ylab = 'TDN (mg N/g soil)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(mb.tdn)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment) %>%
  shapiro_test(mb.tdn)
min(normality$p)

ggqqplot(June2021, "mb.tdn", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")
##the points fall along the line in the qq plot

###ANOVA
JuneTDN.aov <- aov(mb.tdn~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JuneTDN.aov)
##No significance from within effects

res.aov <- anova_test(
  data = June2021, dv = mb.tdn, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance

##June TDN plot
June_TDN <- ggline(June2021, x="date", y = "mb.tdn", color = "treatment", title = "",
                   add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                   size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  labs(x = '', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =30, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=30, label = "Post", color = 'black', size = 9) +
  ylim(0,30)


##Todd's Meadow
summary <- TM2021 %>%
  group_by(sampling.date) %>%
  get_summary_stats(mb.tdn, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "mb.tdn",
          color = "treatment", palette = 'blue', shape = "treatment", add = c("mean_se", "jitter"),
          title = "TDN by Treatment over June 2021 Pulse Experiment", ylab = 'TDN (mg N/g soil)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)

outliers <- TM2021 %>%
  group_by(sampling.date) %>%
  identify_outliers(mb.tdn)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(sampling.date) %>%
  shapiro_test(mb.tdn)
min(normality$p)

ggqqplot(TM2021, "mb.tdn", ggtheme = theme_bw()) 
##normality check p>0.05; qqplot shows majority of points are along the regression line; accept normality check 

###ANOVA
TM_TDN.aov <- aov(mb.tdn~date + Error(core.id/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TM_TDN.aov)
##No significance from within effects

res.aov <- anova_test(
  data = TM2021, dv = mb.tdn, wid = core.id,
  between = c(date)
)
get_anova_table(res.aov)
##No significance

##June TDN plot
TM_TDN <- ggline(TM2021, x="date", y = "mb.tdn", color = "treatment", title = "",
                 add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                 size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = 'blue') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  labs(x = '', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =30, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=30, label = "Post", color = 'black', size = 9) +
  ylim(0,30)


##Between Sampling Campaigns
summary <- all_pulse %>%
  group_by(sampling.campaign) %>%
  get_summary_stats(mb.tdn, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "mb.tdn",
          color = "sampling.campaign", palette = c('grey', 'orange', 'blue'), shape = "sampling.campaign", add = c("mean_se", "jitter"),
          title = "DOC between Experimental and Natural Pulses", ylab = 'DOC (mg C/g soil)', 
          xlab = "", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(sampling.campaign) %>%
  identify_outliers(mb.tdn)
outliers
##extreme outliers detected but no variable removed due to sample size

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(sampling.campaign) %>%
  shapiro_test(mb.tdn)
min(normality$p)

ggqqplot(all_pulse, "mb.tdn", ggtheme = theme_bw()) +
  facet_grid(~sampling.campaign, labeller = "label_both")
##majority of points fall along the line

###ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = mb.tdn, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)
##no significance


TDN_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "mb.tdn",
                                   color = "black", add = c("mean_se", "jitter"),
                                   fill = 'sampling.campaign',
                                   palette =c('grey', 'orange', 'blue'),
                                   legend.title = "Sampling Campaign", legend = "bottom",
                                   title = "Microbial Biomass N Between Sampling Campaigns", ylab = 'Microbial Biomass N (mg N/g soil)', 
                                   xlab = "", bxp.errorbar = TRUE)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 26, face = 'bold', color = 'black'),
        legend.key.size = unit(1, 'cm'), legend.text = element_text(size = 20)) +
  ylim(0,30)




# Microbial Biomass C:N ---------------------------------------------------
##Sept 2020
summary <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(CN, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "CN",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "C:N by Treatment over Sept 2020 Pulse Experiment", ylab = 'C:N', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(CN)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(CN)
min(normality$p)

ggqqplot(Sept2020, "CN", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")
##the points fall along the line in the qq plot

###ANOVA
SeptCN.aov <- aov(CN~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SeptCN.aov)
##Treatment Effect and treatment:date interaction effect

##Treatment post-hoc
emm <- emmeans(SeptCN.aov, ~treatment|date)
p <- pairs(emm)

##interaction post hoc
##post hoc interaction terms
emm <- emmeans(SeptCN.aov, ~date|treatment)
p <- pairs(emm)
##significance shown within sampling groups at different days


res.aov <- anova_test(
  data = Sept2020, dv = CN, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significanc on between terms


Sept_CN <- ggline(Sept2020, x="date", y = "CN", color = "treatment", title = "",
                  add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                  size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  labs(x = 'Day of Pulse', y= 'C:N') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =80, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=80, label = "Post", color = 'black', size = 9) +
  ylim(0,80)


##June 2021
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(CN, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "CN",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "C:N by Treatment over June 2021 Pulse Experiment", ylab = 'C:N', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(CN)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(CN)
min(normality$p)

ggqqplot(June2021, "CN", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")
#majority of points fall on line

###ANOVA
JuneCN.aov <- aov(CN~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JuneCN.aov)
##treatment:date interaction effect

##interaction post hoc
##post hoc interaction terms
emmip(JuneCN.aov, treatment ~ date) ##visualize data
emm <- emmeans(JuneCN.aov, ~date|treatment)
p <- pairs(emm)

emm <- emmeans(JuneCN.aov, ~treatment|date)
p <- pairs(emm)
##significance shown within sampling groups at different days

res.aov <- anova_test(
  data = June2021, dv = CN, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance


June_CN <- ggline(June2021, x="date", y = "CN", color = "treatment", title = "",
                  add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                  size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  labs(x = 'Day of Pulse', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =80, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=80, label = "Post", color = 'black', size = 9) +
  ylim(0, 80)


##Todd's Meadow
summary <- TM2021_1 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(CN, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021_1, x = "date", y = "CN",
          color = "treatment", palette = 'blue', shape = "treatment", add = c("mean_se", "jitter"),
          title = "Carbon:Nitrogen by Treatment over Todd's Meadow Natural Pulse", ylab = 'C:N', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021_1 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(CN)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021_1 %>%
  shapiro_test(CN)
min(normality$p)

ggqqplot(TM2021_1, "CN", ggtheme = theme_bw()) 


###ANOVA
TMCN.aov <- aov(CN~date + Error(core.id/(date)), data = TM2021_1) ##check interaction terms for sampling date and treatment
summary(TMCN.aov)
##date effect

##interaction post hoc
##post hoc interaction terms
emmip(TMCN.aov,  ~ date) ##visualize data
emm <- emmeans(TMCN.aov, ~date)
p <- pairs(emm)
##significance shown within sampling groups at different days

res.aov <- anova_test(
  data = TM2021_1, dv = CN, wid = core.id,
  between = date
)
get_anova_table(res.aov)
##significane of treatment term


TM_CN <- ggline(TM2021_1, x="date", y = "CN", color = "treatment", title = "",
                add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = 'blue') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  labs(x = 'Day of Pulse', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =80, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=80, label = "Post", color = 'black', size = 9) +
  ylim(0,80) +
  geom_bracket(
    xmin = c("0","0"), xmax = c("1", "3"),
    y.position = c(61, 69),
    label = c("*","*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)


##Between Sampling Campaigns
summary <- all_pulse1 %>%
  group_by(sampling.campaign) %>%
  get_summary_stats(CN, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse1, x = "sampling.campaign", y = "CN",
          color = "sampling.campaign", palette = c('grey', 'orange', 'blue'), shape = "sampling.campaign", add = c("mean_se", "jitter"),
          title = "CN between Experimental and Natural Pulse", ylab = 'C:N', 
          xlab = "", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse1 %>%
  group_by(sampling.campaign) %>%
  identify_outliers(CN)
outliers
##extreme outliers detected but no variable removed due to sample size

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse1 %>%
  group_by(sampling.campaign) %>%
  shapiro_test(CN)
min(normality$p)

ggqqplot(all_pulse1, "CN", ggtheme = theme_bw()) +
  facet_grid(~sampling.campaign, labeller = "label_both")
##majority of points fall along the line

###ANOVA
res.aov <- anova_test(
  data = all_pulse1, dv = CN, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)
##no sampling campaign signifcance


CN_sampling.campaign <- ggboxplot(all_pulse1, x = "sampling.campaign", y = "CN",
                                  color = "black", add = c("mean_se", "jitter"),
                                  fill = 'sampling.campaign',
                                  palette =c('grey', 'orange', 'blue'),
                                  legend.title = "Sampling Campaign", legend = "bottom",
                                  title = "Microbial Biomass C:N Ratio Between Sampling Campaigns", ylab = 'C:N', 
                                  xlab = "", bxp.errorbar = TRUE)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 26, face = 'bold', color = 'black'),
        legend.key.size = unit(1, 'cm'), legend.text = element_text(size = 20))+
  ylim(0,75)


ggarrange(Sept_DOC, June_DOC, TM_DOC,
          Sept_TDN, June_TDN, TM_TDN,
          Sept_CN, June_CN, TM_CN,
          nrow = 3, ncol = 3, legend = "none", labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'), 
          font.label = list(size = 30))

ggarrange(DOC_sampling.campaign, TDN_sampling.campaign, CN_sampling.campaign,
          nrow = 1, ncol = 3, legend = "bottom", common.legend = TRUE, labels = c('A', 'B', 'C'), 
          font.label = list(size = 30))

