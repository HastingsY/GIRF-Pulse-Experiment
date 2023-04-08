##Organic Matter, Protein Concentration, N pools and fluxes
##Yvette Hastings
##May 24, 2022
##Updated: April 2, 2023

## Data analysis sections; quick access using shift + alt + J
## 1. Data clean-up
## 2. % organic matter
## 3. Native Protein Concentration
## 4. Organic N
## 5. Inorganic N
## 6. Protein Turnover
## 7. Proteolytic Rate
## 8. Final figures

library(dplyr)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(stringr)
library(rstatix)
library(ggpubr)
library(lubridate)
library(corrplot)
library(reshape)
library(plyr)
library(datarium)
library(sjPlot)
library(ggpmisc)
library(agricolae)
library(emmeans)
library(smatr)
library(jtools)


# Data clean-up -----------------------------------------------------------

##load data
all_pulse <- read_excel("Data/2020-2021 Soils Master File.xlsx", sheet = 'Pulse Data Masterfile (R)') 

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

June2021 <- filter(all_pulse, sampling.campaign == "June 2021 GIRF Pulse")  
Sept2020 <- filter(all_pulse, sampling.campaign == "GIRF September 2020 Pulse")
TM2021 <- filter(all_pulse, sampling.campaign == "June 2021 TM Natural Pulse")



# % organic matter --------------------------------------------------------


##Sept 2020

summary <- Sept2020 %>%
  group_by(treatment, date) %>%
  get_summary_stats(percent.organic.matter, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "percent.organic.matter",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "%OM by Treatment over Sept 2020 Pulse Experiment", ylab = '% of Dry Mass', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(percent.organic.matter)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(percent.organic.matter)
min(normality$p)

ggqqplot(Sept2020, "percent.organic.matter", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
SeptOM.aov <- aov(percent.organic.matter~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SeptOM.aov)
##No significance from within effects


res.aov <- anova_test(
  data = Sept2020, dv = percent.organic.matter, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##Between effect significance by date and treatment

##post hoc on date effect
Sept2020 %>%
  pairwise_t_test(
    percent.organic.matter ~ date,  
    p.adjust.method = "bonferroni"
  )

##post hoc on treatment effect
Sept2020 %>%
  pairwise_t_test(
    percent.organic.matter ~ treatment,  
    p.adjust.method = "bonferroni"
  )


Sept_OM <- ggline(Sept2020, x="date", y = "percent.organic.matter", color = "treatment", title = "September 2020",
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
  labs(x = '', y= '% of Dry Mass') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =15, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=15, label = "Post", color = 'black', size = 9)+ylim (0,15) +
  geom_bracket(
    xmin = c("5"), xmax = c("16"),
    y.position = c(11),
    label = c("*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)

Sept_OM_treatment <- ggboxplot(Sept2020, x = "treatment", y = "percent.organic.matter",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = "September 2020", ylab = '% of Dry Mass', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  geom_bracket(
    xmin = c("Diverse","Diverse"), xmax = c("Grass", "Mulch"),
    y.position = c(11, 9),
    label = c("*","*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,15)

##June 2021
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(percent.organic.matter, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "percent.organic.matter",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "%OM by Treatment over June 2021 Pulse Experiment", ylab = '% of Dry Mass', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(percent.organic.matter)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(percent.organic.matter)
min(normality$p)

ggqqplot(Sept2020, "percent.organic.matter", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
JuneOM.aov <- aov(percent.organic.matter~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JuneOM.aov)
##No significance from within effects


res.aov <- anova_test(
  data = June2021, dv = percent.organic.matter, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance

June_OM <- ggline(June2021, x="date", y = "percent.organic.matter", color = "treatment", title = "June 2021",
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  labs(x = '', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =15, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.6, y=15, label = "Post", color = 'black', size = 9)+ylim (0,15) 

June_OM_treatment <- ggboxplot(June2021, x = "treatment", y = "percent.organic.matter",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = "June 2021", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  ylim(0,15)

##Todd's Meadow
summary <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(percent.organic.matter, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "percent.organic.matter",
          color = "treatment", palette = c('blue'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "June 2021 Experimental Pulse % Organic Matter", ylab = '% of Dry Mass', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(percent.organic.matter)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(percent.organic.matter)
min(normality$p)

ggqqplot(TM2021, "percent.organic.matter", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
TMOM.aov <- aov(percent.organic.matter~date + Error(plot.number/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TMOM.aov)
##No significance from within effects


res.aov <- anova_test(
  data = TM2021, dv = percent.organic.matter, wid = plot.number,
  between = c(date)
)
get_anova_table(res.aov)
##No significance

TM_OM <- ggline(TM2021, x="date", y = "percent.organic.matter", color = "treatment", title = "Todd's Meadow",
                add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('blue')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  labs(x = '', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =15, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=15, label = "Post", color = 'black', size = 9) +
  ylim(0,15)

TM_OM_treatment <- ggboxplot(TM2021, x = "treatment", y = "percent.organic.matter",
                             color = "black", add = c("mean_se", "jitter"),
                             fill =c('blue'),
                             title = "Todd's Meadow", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  ylim(0,15)


##Between Sampling Campaigns
summary <- all_pulse %>%
  group_by(sampling.campaign) %>%
  get_summary_stats(percent.organic.matter, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "percent.organic.matter",
          color = "treatment", shape = "treatment", add = c("mean_se", "jitter"),
          palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue'),
          title = "%OM Between Sampling Campaigns", ylab = '% of Dry Mass', 
          xlab = "Sampling Campaign", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  identify_outliers(percent.organic.matter)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  shapiro_test(percent.organic.matter)
min(normality$p)


##ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = percent.organic.matter, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)


all_pulse %>%
  pairwise_t_test(
    percent.organic.matter ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )

##plot results
OM_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "percent.organic.matter",
                                  color = "black", add = c("mean_se", "jitter"),
                                  fill =c('grey', 'orange', 'blue'),
                                  title = "Between Sampling Campaigns", ylab = '', 
                                  xlab = "", bxp.errorbar = TRUE)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  scale_x_discrete(labels=c("GIRF September 2020 Pulse" = "GIRF Sept 2020", "June 2021 GIRF Pulse" = "GIRF June 2021", "June 2021 TM Natural Pulse" = "Reference June 2021")) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse", "GIRF September 2020 Pulse", "June 2021 GIRF Pulse"), xmax = c("June 2021 GIRF Pulse", "June 2021 TM Natural Pulse", "June 2021 TM Natural Pulse"), 
    y.position = c(14.5, 13, 11.5),
    label = c("**", "***", "***"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  ylim(0,15)



# Native Protein Concentration --------------------------------------------


##Sept 2020
summary <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(Native.Protein.Concentration, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "Native.Protein.Concentration",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = 'September 2020 Experimental Pulse Native Protein', ylab = 'ug protein/g soil', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(Native.Protein.Concentration)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(Native.Protein.Concentration)
min(normality$p)

ggqqplot(Sept2020, "Native.Protein.Concentration", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
SeptNP.aov <- aov(Native.Protein.Concentration~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SeptNP.aov)
##No significance


res.aov <- anova_test(
  data = Sept2020, dv = Native.Protein.Concentration, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##Date and treatment effect

Sept2020 %>%
  pairwise_t_test(
    Native.Protein.Concentration ~ date,  
    p.adjust.method = "bonferroni"
  )

Sept2020 %>%
  pairwise_t_test(
    Native.Protein.Concentration ~ treatment,  
    p.adjust.method = "bonferroni"
  )


Sept_NP <- ggline(Sept2020, x="date",  y = "Native.Protein.Concentration", color = "treatment", title = '',
                  add = c("mean_se"), palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                  legend.title = "Treatment", legend = 'none', size =1, shape = "treatment", point.size = 5)+
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =40, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=40, label = "Post", color = 'black', size = 9) +
  ylab('ug protein/g soil')+
  xlab("Day of Pulse") +
  ylim(0,40)+
  geom_bracket(
    xmin = c("3"), xmax = c("16"),
    y.position = c(30),
    label = c("**"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)

Sept_NP_treatment <- ggboxplot(Sept2020, x = "treatment", y = "Native.Protein.Concentration",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = '', ylab = 'ug protein/g soil', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  geom_bracket(
    xmin = c("Diverse"), xmax = c("Grass"),
    y.position = c(28),
    label = c("*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  ylim(0, 40)

##June 2021
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(Native.Protein.Concentration, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "Native.Protein.Concentration",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "June 2021 Native Protein", ylab = 'ug Protein/g soil', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(Native.Protein.Concentration)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(Native.Protein.Concentration)
min(normality$p)

ggqqplot(June2021, "Native.Protein.Concentration", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
JuneNP.aov <- aov(Native.Protein.Concentration~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JuneNP.aov)
##No significance from within effects


res.aov <- anova_test(
  data = June2021, dv = Native.Protein.Concentration, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##date effect

##post hoc on date effect
June2021 %>%
  pairwise_t_test(
    Native.Protein.Concentration ~ treatment,   
    p.adjust.method = "bonferroni"
  )


June_NP <- ggline(June2021, x="date",  y = "Native.Protein.Concentration", color = "treatment", title = '',
                  add = c("mean_se"), palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                  legend.title = "Treatment", legend = 'none', size =1, shape = "treatment", point.size = 5)+
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =40, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=40, label = "Post", color = 'black', size = 9) +
  ylab('')+
  xlab("Day of Pulse") +
  ylim(0,40)


June_NP_treatment <- ggboxplot(June2021, x = "treatment", y = "Native.Protein.Concentration",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = "", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  geom_bracket(
    xmin = c("Grass"), xmax = c("Mulch"),
    y.position = c(32),
    label = c("***"),
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  ylim(0, 40)

##Todd's Meadow
summary <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(Native.Protein.Concentration, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "Native.Protein.Concentration",
          color = "treatment", palette = c('blue'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "2021 Todd's Meadow Native Protein", ylab = 'ug protein/g soil', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(Native.Protein.Concentration)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(Native.Protein.Concentration)
min(normality$p)

ggqqplot(TM2021, "Native.Protein.Concentration", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
TMNP.aov <- aov(Native.Protein.Concentration~date + Error(plot.number/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TMNP.aov)
##No significance from within effects


res.aov <- anova_test(
  data = TM2021, dv = Native.Protein.Concentration, wid = plot.number,
  between = c(date)
)
get_anova_table(res.aov)
##Date effect

##post hoc on treatment effect
TM2021 %>%
  pairwise_t_test(
    Native.Protein.Concentration ~ date,  
    p.adjust.method = "bonferroni"
  )

TM_NP <- ggline(TM2021, x="date",  y = "Native.Protein.Concentration", color = "treatment", title = "",
                add = c("mean_se"), palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                legend.title = "Treatment", legend = 'none', size =1, shape = "treatment", point.size = 5)+
  scale_color_manual(values = c('blue')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =40, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=40, label = "Post", color = 'black', size = 9) +
  ylab('')+
  xlab("Day of Pulse") +
  geom_bracket(
    xmin = c("0"), xmax = c("3"),
    y.position = c(36),
    label = c("*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  ylim(0, 40)

TM_NP_treatment <- ggboxplot(TM2021, x = "treatment", y = "Native.Protein.Concentration",
                             color = "black", add = c("mean_se", "jitter"),
                             fill =c('blue'),
                             title = "Todd's Meadow", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))


##Between sampling campaigns
summary <- all_pulse %>%
  group_by(sampling.campaign) %>%
  get_summary_stats(Native.Protein.Concentration, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "Native.Protein.Concentration",
          color = "treatment", shape = "treatment", add = c("mean_se", "jitter"),
          palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue'),
          title = "Native Protein Between Sampling Campaigns", ylab = 'ug protein/g soil', 
          xlab = "Sampling Campaign", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  identify_outliers(Native.Protein.Concentration)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(sampling.campaign) %>%
  shapiro_test(Native.Protein.Concentration)
min(normality$p)


##ANOVA
allNP.aov <- aov(Native.Protein.Concentration~sampling.campaign + Error(plot.number/(sampling.campaign)), data = all_pulse) ##check interaction terms for sampling date and treatment
summary(allNP.aov)


res.aov <- anova_test(
  data = all_pulse, dv = Native.Protein.Concentration, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)


all_pulse %>%
  pairwise_t_test(
    Native.Protein.Concentration ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )


##plot results
NP_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "Native.Protein.Concentration",
                                  color = "black", add = c("mean_se", "jitter"),
                                  fill =c('grey', 'orange', 'blue'),
                                  title = "", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  ylim(0, 45) +
  scale_x_discrete(labels=c("GIRF September 2020 Pulse" = "GIRF Sept 2020", "June 2021 GIRF Pulse" = "GIRF June 2021", "June 2021 TM Natural Pulse" = "Reference June 2021")) +
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse", "June 2021 GIRF Pulse"), xmax = c("June 2021 TM Natural Pulse", "June 2021 TM Natural Pulse"), 
    y.position = c(44,40),
    label = c("***", "***"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5)



# Organic N ---------------------------------------------------------------



##Sept 2020
summary <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(organic.nitrogen, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "organic.nitrogen",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "Organic Nitrogen by Treatment over Sept 2020 Pulse Experiment", ylab = 'mg N/g soil', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(organic.nitrogen)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(organic.nitrogen)
min(normality$p)

ggqqplot(Sept2020, "organic.nitrogen", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
SeptON.aov <- aov(organic.nitrogen~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SeptON.aov)
##No significance from within effects


res.aov <- anova_test(
  data = Sept2020, dv = organic.nitrogen, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##Between effect significance by date and treatment

##post hoc on date effect
Sept2020 %>%
  pairwise_t_test(
    organic.nitrogen ~ date,  
    p.adjust.method = "bonferroni"
  )


##post hoc on treatment effect
Sept2020 %>%
  pairwise_t_test(
    organic.nitrogen ~ treatment,  
    p.adjust.method = "bonferroni"
  )


Sept_ON <- ggline(Sept2020, x="date", y = "organic.nitrogen", color = "treatment", title = "September 2020",
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  labs(x = '', y= 'mg N/g soil') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =5, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=5, label = "Post", color = 'black', size = 9)+ylim (0,5) +
  geom_bracket(
    xmin = c("3", "5", "9"), xmax = c("16", "16", "16"),
    y.position = c(3, 2.5, 2),
    label = c("*", "**", "**"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)

Sept_ON_treatment <- ggboxplot(Sept2020, x = "treatment", y = "organic.nitrogen",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = "September 2020", ylab = 'mg N/g soil', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  geom_bracket(
    xmin = c("Diverse"), xmax = c("Grass"),
    y.position = c(2.5),
    label = c("*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,5)

##June 2021
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(organic.nitrogen, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "organic.nitrogen",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "Organic Nitrogen by Treatment over June 2021 Pulse Experiment", ylab = '% of Dry Mass', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(organic.nitrogen)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(organic.nitrogen)
min(normality$p)

ggqqplot(June2021, "organic.nitrogen", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
JuneON.aov <- aov(organic.nitrogen~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JuneON.aov)
##No significance from within effects


res.aov <- anova_test(
  data = June2021, dv = organic.nitrogen, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance

June_ON <- ggline(June2021, x="date", y = "organic.nitrogen", color = "treatment", title = "June 2021",
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  labs(x = '', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =5, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=5, label = "Post", color = 'black', size = 9)+ylim (0,5) 

June_ON_treatment <- ggboxplot(June2021, x = "treatment", y = "organic.nitrogen",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = "June 2021", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

##Todd's Meadow
summary <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(organic.nitrogen, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "organic.nitrogen",
          color = "treatment", palette = c('blue'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "Organic Nitrogen by Treatment over 2021 Todd's Meadow Natural Pulse", ylab = 'mg N/g soil', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(organic.nitrogen)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(organic.nitrogen)
min(normality$p)

ggqqplot(TM2021, "organic.nitrogen", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
TMON.aov <- aov(organic.nitrogen~date + Error(plot.number/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TMON.aov)
##No significance from within effects


res.aov <- anova_test(
  data = TM2021, dv = organic.nitrogen, wid = plot.number,
  between = c(date)
)
get_anova_table(res.aov)
##No significance

TM_ON <- ggline(TM2021, x="date", y = "organic.nitrogen", color = "treatment", title = "Todd's Meadow",
                add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('blue')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  labs(x = '', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =5, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=5, label = "Post", color = 'black', size = 9) + ylim(0, 5)

TM_ON_treatment <- ggboxplot(TM2021, x = "treatment", y = "inorganic.nitrogen",
                             color = "black", add = c("mean_se", "jitter"),
                             fill =c('blue'),
                             title = "Todd's Meadow", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))


##Between Sampling Campaigns
summary <- all_pulse %>%
  group_by(sampling.campaign) %>%
  get_summary_stats(organic.nitrogen, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "organic.nitrogen",
          color = "treatment", shape = "treatment", add = c("mean_se", "jitter"),
          palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue'),
          title = "ON Between Sampling Campaigns", ylab = 'mg N/g soil', 
          xlab = "Sampling Campaign", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  identify_outliers(organic.nitrogen)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(sampling.campaign) %>%
  shapiro_test(organic.nitrogen)
min(normality$p)


##ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = organic.nitrogen, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)


all_pulse %>%
  pairwise_t_test(
    organic.nitrogen ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )


##plot results
ON_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "organic.nitrogen",
                                  color = "black", add = c("mean_se", "jitter"),
                                  fill =c('grey', 'orange', 'blue'),
                                  title = "Between Sampling Campaigns", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  ylim(0, 5) +
  scale_x_discrete(labels=c("GIRF September 2020 Pulse" = "GIRF Sept 2020", "June 2021 GIRF Pulse" = "GIRF June 2021", "June 2021 TM Natural Pulse" = "Reference June 2021")) +
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse", "June 2021 GIRF Pulse"), xmax = c("June 2021 TM Natural Pulse", "June 2021 TM Natural Pulse"), 
    y.position = c(4.9, 4.4),
    label = c("***", "***"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5)




# Inorganic N -------------------------------------------------------------



##Sept 2020
summary <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(inorganic.nitrogen, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "inorganic.nitrogen",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "Inorganic Nitrogen by Treatment over Sept 2020 Pulse Experiment", ylab = 'mg N/g soil', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(inorganic.nitrogen)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(inorganic.nitrogen)
min(normality$p)

ggqqplot(Sept2020, "inorganic.nitrogen", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
SeptIN.aov <- aov(inorganic.nitrogen~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SeptIN.aov)
##Significan of treatment 


res.aov <- anova_test(
  data = Sept2020, dv = inorganic.nitrogen, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance

emm <- emmeans(SeptIN.aov, ~treatment)
p <- pairs(emm)


Sept_IN <- ggline(Sept2020, x="date", y = "inorganic.nitrogen", color = "treatment", title = "",
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  labs(x = 'Day of Pulse', y= 'mg N/g soil') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =0.2, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=0.2, label = "Post", color = 'black', size = 9)+ylim (0,0.2) 

Sept_IN_treatment <- ggboxplot(Sept2020, x = "treatment", y = "inorganic.nitrogen",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = "", ylab = 'mg N/g soil', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  geom_bracket(
    xmin = c("Diverse"), xmax = c("Grass"),
    y.position = c(0.1),
    label = c("*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,0.2)

##June 2021
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(inorganic.nitrogen, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "inorganic.nitrogen",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "Inorganic Nitrogen by Treatment over June 2021 Pulse Experiment", ylab = 'mg N/g soil', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(inorganic.nitrogen)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(inorganic.nitrogen)
min(normality$p)

ggqqplot(June2021, "inorganic.nitrogen", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
JuneIN.aov <- aov(inorganic.nitrogen~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JuneIN.aov)
##No significance from within effects


res.aov <- anova_test(
  data = June2021, dv = inorganic.nitrogen, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance


June_IN <- ggline(June2021, x="date", y = "inorganic.nitrogen", color = "treatment", title = "",
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  labs(x = 'Day of Pulse', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =0.2, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=0.2, label = "Post", color = 'black', size = 9)+ylim (0,0.2) 

June_IN_treatment <- ggboxplot(June2021, x = "treatment", y = "inorganic.nitrogen",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = "", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

##Todd's Meadow
summary <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(inorganic.nitrogen, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "inorganic.nitrogen",
          color = "treatment", palette = c('blue'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "Inorganic Nitrogen by Treatment over 2021 Todd's Meadow Natural Pulse", ylab = 'mg N/g soil', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(inorganic.nitrogen)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(inorganic.nitrogen)
min(normality$p)

ggqqplot(TM2021, "inorganic.nitrogen", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
TMIN.aov <- aov(inorganic.nitrogen~date + Error(plot.number/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TMIN.aov)
##No significance from within effects


res.aov <- anova_test(
  data = TM2021, dv = inorganic.nitrogen, wid = plot.number,
  between = c(date)
)
get_anova_table(res.aov)
##Date effect

##post hoc on treatment effect
TM2021 %>%
  pairwise_t_test(
    inorganic.nitrogen ~ date,  
    p.adjust.method = "bonferroni"
  )

TM_IN <- ggline(TM2021, x="date", y = "inorganic.nitrogen", color = "treatment", title = "",
                add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('blue'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  labs(x = 'Day of Pulse', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =0.2, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=0.2, label = "Post", color = 'black', size = 9) +
  geom_bracket(
    xmin = c("0", "0"), xmax = c("1", "3"),
    y.position = c(0.06, 0.03),
    label = c("*", "**"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,0.2)

TM_IN_treatment <- ggboxplot(TM2021, x = "treatment", y = "inorganic.nitrogen",
                             color = "black", add = c("mean_se", "jitter"),
                             fill =c('blue'),
                             title = "", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))


##Between Sampling Campaigns
summary <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  get_summary_stats(inorganic.nitrogen, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "inorganic.nitrogen",
          color = "treatment", shape = "treatment", add = c("mean_se", "jitter"),
          palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue'),
          title = "IN Between Sampling Campaigns", ylab = 'mg N/g soil', 
          xlab = "Sampling Campaign", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  identify_outliers(inorganic.nitrogen)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(sampling.campaign) %>%
  shapiro_test(inorganic.nitrogen)
min(normality$p)


##ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = inorganic.nitrogen, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)


##plot results
IN_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "inorganic.nitrogen",
                                  color = "black", add = c("mean_se", "jitter"),
                                  fill =c('grey', 'orange', 'blue'),
                                  title = "", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  scale_x_discrete(labels=c("GIRF September 2020 Pulse" = "GIRF Sept 2020", "June 2021 GIRF Pulse" = "GIRF June 2021", "June 2021 TM Natural Pulse" = "Reference June 2021"))+
  ylim(0, 0.2) 



# Protein Turnover --------------------------------------------------------



##Sept 2020
summary <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(potential.turnover.rate, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "potential.turnover.rate",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = 'September 2020 Experimental Pulse Potential Protein Turnover Rate', ylab = 'Time (hours)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(potential.turnover.rate)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(potential.turnover.rate)
min(normality$p)

ggqqplot(Sept2020, "potential.turnover.rate", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
SeptPT.aov <- aov(potential.turnover.rate~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SeptPT.aov)
##No significance


res.aov <- anova_test(
  data = Sept2020, dv = potential.turnover.rate, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance


Sept_PT <- ggline(Sept2020, x="date",  y = "potential.turnover.rate", color = "treatment", title = 'September 2020',
                  add = c("mean_se"), palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                  legend.title = "Treatment", legend = 'none', size =1, shape = "treatment", point.size = 5)+
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =1, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=1, label = "Post", color = 'black', size = 9) +
  ylab('Time (hours)')+
  xlab("Day of Pulse") +
  ylim(0,1)

Sept_PT_treatment <- ggboxplot(Sept2020, x = "treatment", y = "potential.turnover.rate",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = '', ylab = 'Time (hours)', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

##June 2021
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(potential.turnover.rate, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "potential.turnover.rate",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "June 2021 Potential Protein Turnover Rate", ylab = 'Time (hours)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(potential.turnover.rate)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(potential.turnover.rate)
min(normality$p)

ggqqplot(Sept2020, "potential.turnover.rate", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
JunePT.aov <- aov(potential.turnover.rate~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JunePT.aov)
##No significance from within effects


res.aov <- anova_test(
  data = June2021, dv = potential.turnover.rate, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance


June_PT <- ggline(June2021, x="date",  y = "potential.turnover.rate", color = "treatment", title = 'June 2021',
                  add = c("mean_se"), palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                  legend.title = "Treatment", legend = 'none', size =1, shape = "treatment", point.size = 5)+
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =1, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=1, label = "Post", color = 'black', size = 9) +
  ylab('')+
  xlab("Day of Pulse") +
  ylim(0,1)

June_PT_treatment <- ggboxplot(June2021, x = "treatment", y = "potential.turnover.rate",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = "", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  ylim(0, 1)

##TM
summary <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(potential.turnover.rate, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "potential.turnover.rate",
          color = "treatment", palette = c('blue'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "2021 Todd's Meadow Potential Protein Turnover Rate", ylab = 'Time (hours)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(potential.turnover.rate)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(potential.turnover.rate)
min(normality$p)

ggqqplot(TM2021, "potential.turnover.rate", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
TMPT.aov <- aov(potential.turnover.rate~date + Error(plot.number/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TMPT.aov)
##No significance from within effects


res.aov <- anova_test(
  data = TM2021, dv = potential.turnover.rate, wid = plot.number,
  between = c(date)
)
get_anova_table(res.aov)
##No significance


TM_PT <- ggline(TM2021, x="date",  y = "potential.turnover.rate", color = "treatment", title = "Todd's Meadow",
                add = c("mean_se"), palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                legend.title = "Treatment", legend = 'none', size =1, shape = "treatment", point.size = 5)+
  scale_color_manual(values = c('blue')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =1, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=1, label = "Post", color = 'black', size = 9) +
  ylab('')+
  xlab("Day of Pulse")+
  ylim(0,1)

TM_PT_treatment <- ggboxplot(TM2021, x = "treatment", y = "potential.turnover.rate",
                             color = "black", add = c("mean_se", "jitter"),
                             fill =c('blue'),
                             title = "", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))


##Between Sampling Campaigns
summary <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  get_summary_stats(potential.turnover.rate, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "potential.turnover.rate",
          color = "treatment", shape = "treatment", add = c("mean_se", "jitter"),
          palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue'),
          title = "Potential Protein Turnover Rate Between Sampling Campaigns", ylab = 'Time (hours)', 
          xlab = "Sampling Campaign", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  identify_outliers(potential.turnover.rate)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(sampling.campaign) %>%
  shapiro_test(potential.turnover.rate)
min(normality$p)


##ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = potential.turnover.rate, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)


all_pulse %>%
  pairwise_t_test(
    potential.turnover.rate ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )


##plot results
PT_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "potential.turnover.rate",
                                  color = "black", add = c("mean_se", "jitter"),
                                  fill =c('grey', 'orange', 'blue'),
                                  title = "Between Sampling Campaigns", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  ylim(0, 1) +
  scale_x_discrete(labels=c("GIRF September 2020 Pulse" = "GIRF Sept 2020", "June 2021 GIRF Pulse" = "GIRF June 2021", "June 2021 TM Natural Pulse" = "Reference June 2021")) +
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse"), xmax = c("June 2021 GIRF Pulse"), 
    y.position = c(0.8),
    label = c("**"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5)



# Proteolytic Rate --------------------------------------------------------



##Sept 2020
summary <- all_pulse %>%
  group_by(sampling.campaign, date) %>%
  get_summary_stats(potential.net.proteolytic.rate, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "potential.net.proteolytic.rate",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "Proteoltyic Rate over Sept 2020 Pulse Experiment", ylab = '?g/g soil per hour', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(potential.net.proteolytic.rate)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(potential.net.proteolytic.rate)
min(normality$p)

ggqqplot(Sept2020, "potential.net.proteolytic.rate", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
SeptPR.aov <- aov(potential.net.proteolytic.rate~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SeptPR.aov)
##No significance from within effects


res.aov <- anova_test(
  data = Sept2020, dv = potential.net.proteolytic.rate, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance

Sept_PR <- ggline(Sept2020, x="date", y = "potential.net.proteolytic.rate", color = "treatment", title = "September 2020",
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  labs(x = 'Day of Pulse', y= 'µg/g soil per hour') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=60, label = "Post", color = 'black', size = 9)+ylim (0,60)

Sept_PR_treatment <- ggboxplot(Sept2020, x = "treatment", y = "potential.net.proteolytic.rate",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = "", ylab = 'µg/g soil per hour', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

##June
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(potential.net.proteolytic.rate, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "potential.net.proteolytic.rate",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "Proteoltyic Rate over June 2021 Pulse Experiment", ylab = '?g/g soil per hour', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(potential.net.proteolytic.rate)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(potential.net.proteolytic.rate)
min(normality$p)

ggqqplot(June2021, "potential.net.proteolytic.rate", ggtheme = theme_bw()) +
  facet_grid(~treatment, labeller = "label_both")

###ANOVA
JunePR.aov <- aov(potential.net.proteolytic.rate~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JunePR.aov)
##No significance from within effects


res.aov <- anova_test(
  data = June2021, dv = potential.net.proteolytic.rate, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance

June_PR <- ggline(June2021, x="date", y = "potential.net.proteolytic.rate", color = "treatment", title = "June 2021",
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  labs(x = 'Day of Pulse', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=60, label = "Post", color = 'black', size = 9)+ylim (0,60) 

June_PR_treatment <- ggboxplot(June2021, x = "treatment", y = "potential.net.proteolytic.rate",
                               color = "black", add = c("mean_se", "jitter"),
                               fill =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                               title = "", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

##Todd's Meadow
summary <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(potential.net.proteolytic.rate, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "potential.net.proteolytic.rate",
          color = "treatment", palette = c('blue'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "June 2021 Experimental Pulse Proteolytic Rate", ylab = '?g/g soil per hour', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(potential.net.proteolytic.rate)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(potential.net.proteolytic.rate)
min(normality$p)

ggqqplot(TM2021, "potential.net.proteolytic.rate", ggtheme = theme_bw()) +
  facet_grid(treatment ~ date, labeller = "label_both")

###ANOVA
TMPR.aov <- aov(potential.net.proteolytic.rate~date + Error(plot.number/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TMPR.aov)
##No significance from within effects


res.aov <- anova_test(
  data = TM2021, dv = potential.net.proteolytic.rate, wid = plot.number,
  between = c(date)
)
get_anova_table(res.aov)
##No significance

TM_PR <- ggline(TM2021, x="date", y = "potential.net.proteolytic.rate", color = "treatment", title = "Todd's Meadow",
                add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('blue')) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  labs(x = 'Day of Pulse', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=60, label = "Post", color = 'black', size = 9) +
  ylim(0,60)

TM_PR_treatment <- ggboxplot(TM2021, x = "treatment", y = "potential.net.proteolytic.rate",
                             color = "black", add = c("mean_se", "jitter"),
                             fill =c('blue'),
                             title = "", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  ylim(0,50)


##Between Sampling Campaigns
summary <- all_pulse %>%
  group_by(sampling.campaign) %>%
  get_summary_stats(potential.net.proteolytic.rate, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "potential.net.proteolytic.rate",
          color = "treatment", shape = "treatment", add = c("mean_se", "jitter"),
          palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue'),
          title = "PR Between Sampling Campaigns", ylab = '?g / g soil per hour', 
          xlab = "Sampling Campaign", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  identify_outliers(potential.net.proteolytic.rate)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(sampling.campaign) %>%
  shapiro_test(potential.net.proteolytic.rate)
min(normality$p)


##ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = potential.net.proteolytic.rate, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)


all_pulse %>%
  pairwise_t_test(
    potential.net.proteolytic.rate ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )


##plot results
PR_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "potential.net.proteolytic.rate",
                                  color = "black", add = c("mean_se", "jitter"),
                                  fill =c('grey', 'orange', 'blue'),
                                  title = "Between Sampling Campaigns", ylab = '', 
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
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))+
  scale_x_discrete(labels=c("GIRF September 2020 Pulse" = "GIRF Sept 2020", "June 2021 GIRF Pulse" = "GIRF June 2021", "June 2021 TM Natural Pulse" = "Reference June 2021")) +
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse"), xmax = c("June 2021 GIRF Pulse"), 
    y.position = c(59),
    label = c("*"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  ylim(0,60)



# final figures -----------------------------------------------------------

ggarrange(Sept_OM, June_OM, TM_OM,
          Sept_NP, June_NP, TM_NP,
          nrow = 2, ncol = 3, legend = "none", labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
          font.label = list(size = 30))

ggarrange(Sept_ON, June_ON, TM_ON,
          Sept_IN, June_IN, TM_IN,
          nrow = 2, ncol = 3, legend = "none", labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
          font.label = list(size = 30))

ggarrange(Sept_OM_treatment, June_OM_treatment, OM_sampling.campaign,
          Sept_NP_treatment, June_NP_treatment, NP_sampling.campaign,
          nrow = 2, ncol = 3, legend = "none", labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
          font.label = list(size = 30))

ggarrange(Sept_ON_treatment, June_ON_treatment, ON_sampling.campaign,
          Sept_IN_treatment, June_IN_treatment, IN_sampling.campaign,
          nrow = 2, ncol = 3, legend = "none", labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
          font.label = list(size = 30))

ggarrange(Sept_PR, June_PR, TM_PR,
          Sept_PR_treatment, June_PR_treatment, PR_sampling.campaign,
          nrow = 2, ncol = 3, legend = "none", labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
          font.label = list(size = 30))

ggarrange(Sept_PT, June_PT, TM_PT,
          Sept_PT_treatment, June_PT_treatment, PT_sampling.campaign,
          nrow = 2, ncol = 3, legend = "none", labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
          font.label = list(size = 30))

