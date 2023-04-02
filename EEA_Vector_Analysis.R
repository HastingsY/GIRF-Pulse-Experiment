##Ecoenzyme SMA and Vector Analysis
##Yvette Hastings
##April 17, 2022
##Updated 4/2/2023

## Data analysis sections; quick access using shift + alt + J
## 1. Data clean-up
## 2. Ecoenzyme RMANOVA
## 3. SMA for EEA
## 4. Vector Analysis

##load libraries
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
library(smatr)
library(emmeans)
library(pracma) ##for rad2deg to calculate degree in vector analysis
library(ggtext) ##add R2 syntax to SMA plots


# Data clean-up -----------------------------------------------------------
##load and data formatting
all_pulse <- read_excel("Data/2020-2021 Soils Master File.xlsx", sheet = 'Pulse Data Masterfile (R)') 

all_pulse$lnAP <- log(all_pulse$enzyme.ap)
all_pulse$lnLAP <- log(all_pulse$enzyme.lap)
all_pulse$lnBG <- log(all_pulse$enzyme.bg)
all_pulse$lnPOX <- log(all_pulse$enzyme.pox)
all_pulse$lnBG_lnAP <- all_pulse$lnBG/all_pulse$lnAP
all_pulse$lnBG_lnLAP <- all_pulse$lnBG/all_pulse$lnLAP
all_pulse$length <-sqrt((all_pulse$lnBG_lnAP^2)+(all_pulse$lnBG_lnLAP^2))
all_pulse$angle <- rad2deg(atan2(all_pulse$lnBG_lnAP, all_pulse$lnBG_lnLAP))

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


# all_pulse <- filter(all_pulse, sampling.campaign == "GIRF September 2020 Pulse" & sampling.campign == "June 2021 GIRF Pulse")
June2021 <- filter(all_pulse, sampling.campaign == "June 2021 GIRF Pulse") 
Sept2020 <- filter(all_pulse, sampling.campaign == "GIRF September 2020 Pulse")
TM2021 <- filter(all_pulse, sampling.campaign == "June 2021 TM Natural Pulse")
June2021_1 <- filter(all_pulse, sampling.campaign == "June 2021 GIRF Pulse" | sampling.campaign == "June 2021 TM Natural Pulse")



# Ecoenzyme RMANOVA -------------------------------------------------------
##Sept 2020

##LAP
summary <- Sept2020 %>%
  group_by(sampling.date) %>%
  get_summary_stats(enzyme.lap, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "enzyme.lap",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "LAP activity Sept 2020", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(enzyme.lap)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(enzyme.lap)
min(normality$p)

ggqqplot(Sept2020, "enzyme.lap", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
SetpLAP.aov <- aov(enzyme.lap~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SetpLAP.aov)
##No significance from within effects

res.aov <- anova_test(
  data = Sept2020, dv = enzyme.lap, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##* significance for date

##post hoc on date effect
Sept2020 %>%
  pairwise_t_test(
    enzyme.lap ~ date, 
    p.adjust.method = "bonferroni"
  )


LAP_Sept_plot<- ggline(Sept2020, x="date", y = "enzyme.lap", color = "treatment", title = "September 2020",
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
  labs(x = '', y= expression(LAP~(nmol~g^-1~h^-1))) +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=60, label = "Post", color = 'black', size = 9)+
  geom_bracket(
    xmin = c("0","0", "0", "0"), xmax = c("16", "9", "3", "5"), 
    y.position = c(36,32,28,24),
    label = c("*","***", "*", "***"), #based on p-adjusted 
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,60)

##AP
summary <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(enzyme.ap, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "enzyme.ap",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "AP activity Sept 2020", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(enzyme.ap)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(enzyme.ap)
min(normality$p)

ggqqplot(Sept2020, "enzyme.ap", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
SetpAP.aov <- aov(enzyme.ap~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SetpAP.aov)
##no significance 


res.aov <- anova_test(
  data = Sept2020, dv = enzyme.ap, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##* significance for date

##post hoc on date effect
Sept2020 %>%
  pairwise_t_test(
    enzyme.ap ~ date, 
    p.adjust.method = "bonferroni"
  )

AP_Sept_plot<- ggline(Sept2020, x="date", y = "enzyme.ap", color = "treatment", title = "",
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
  labs(x = 'Day of Pulse', y= expression(AP~(nmol~g^-1~h^-1))) +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=60, label = "Post", color = 'black', size = 9)+
  ylim (0,60)

##BG
summary <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(enzyme.bg, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "enzyme.bg",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "BG activity Sept 2020", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(enzyme.bg)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(enzyme.bg)
min(normality$p)

ggqqplot(Sept2020, "enzyme.bg", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
SetpBG.aov <- aov(enzyme.bg~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SetpBG.aov)
##No significance from within effects

res.aov <- anova_test(
  data = Sept2020, dv = enzyme.bg, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##* significance for date

##post hoc on date effect
Sept2020 %>%
  pairwise_t_test(
    enzyme.bg ~ date, 
    p.adjust.method = "bonferroni"
  )

BG_Sept_plot<- ggline(Sept2020, x="date", y = "enzyme.bg", color = "treatment", title = "September 2020",
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
  labs(x = '', y= expression(BG~(nmol~g^-1~h^-1))) +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=60, label = "Post", color = 'black', size = 9)+
  geom_bracket(
    xmin = c("0", "0", "0"), xmax = c("9", "5", "3"), 
    y.position = c(28,24 , 20),
    label = c("*","*", "*"), #based on p-adjusted
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,60)

##POX
summary <- Sept2020 %>%
  group_by(sampling.date) %>%
  get_summary_stats(enzyme.pox, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "enzyme.pox",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "POX activity Sept 2020", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(enzyme.pox)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(enzyme.pox)
min(normality$p)

ggqqplot(Sept2020, "enzyme.pox", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
SetpPOX.aov <- aov(enzyme.pox~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SetpPOX.aov)
##No significance from within effects

res.aov <- anova_test(
  data = Sept2020, dv = enzyme.pox, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##* significance for date


POX_Sept_plot<- ggline(Sept2020, x="date", y = "enzyme.pox", color = "treatment", title = "",
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
  labs(x = 'Day of Pulse', y= expression(POX~(nmol~g^-1~h^-1))) +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =30000, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=30000, label = "Post", color = 'black', size = 9) +
  ylim(0, 30000)



##June 2021
##LAP
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(enzyme.lap, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "enzyme.lap",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "LAP activity June 2021", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(enzyme.lap)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(enzyme.lap)
min(normality$p)

ggqqplot(June2021, "enzyme.lap", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
JuneLAP.aov <- aov(enzyme.lap~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JuneLAP.aov)
##No significance from within effects

res.aov <- anova_test(
  data = June2021, dv = enzyme.lap, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##* significance for date and treatment

##post hoc on date effect
June2021 %>%
  pairwise_t_test(
    enzyme.lap ~ date, 
    p.adjust.method = "bonferroni"
  )


LAP_June_plot<- ggline(June2021, x="date", y = "enzyme.lap", color = "treatment", title = "June 2021",
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
        legend.title = element_text(size = 16, face = 'bold', color = 'black'),
        legend.key.size = unit(1, 'cm'), legend.text = element_text(size = 20))+
  labs(x = '', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=60, label = "Post", color = 'black', size = 9)+  ##no ns after post hoc
  ylim (0,60)


##June AP
summary <- June2021 %>%
  group_by(sampling.date) %>%
  get_summary_stats(enzyme.ap, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "enzyme.ap",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "AP activity June 2021", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(enzyme.ap)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(enzyme.ap)
min(normality$p)

ggqqplot(June2021, "enzyme.ap", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
JuneAP.aov <- aov(enzyme.ap~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JuneAP.aov)
##No significance from within effects

res.aov <- anova_test(
  data = June2021, dv = enzyme.ap, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##* significance for date and treatment

##post hoc on date effect
p <- June2021 %>%
  pairwise_t_test(
    enzyme.ap ~ date, 
    p.adjust.method = "bonferroni"
  )
##** adj.p.sig 0-3 day


AP_June_plot<- ggline(June2021, x="date", y = "enzyme.ap", color = "treatment", title = "",
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
  geom_vline(xintercept = 1.65, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=60, label = "Post", color = 'black', size = 9)+
  geom_bracket(
    xmin = c("0"), xmax = c("1"), 
    y.position = c(50),
    label = c("**"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,60)

##June BG
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(enzyme.bg, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "enzyme.bg",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "BG activity June 2021", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(enzyme.bg)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(enzyme.bg)
min(normality$p)

ggqqplot(June2021, "enzyme.bg", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
JuneBG.aov <- aov(enzyme.bg~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JuneBG.aov)
##* significance on date

pwc <- June2021 %>%
  pairwise_t_test(
    enzyme.bg ~date, 
    p.adjust.method = 'bonferroni'
  )
pwc

res.aov <- anova_test(
  data = June2021, dv = enzyme.bg, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##No significance


BG_June_plot<- ggline(June2021, x="date", y = "enzyme.bg", color = "treatment", title = "June 2021",
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
  geom_vline(xintercept = 1.65, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=60, label = "Post", color = 'black', size = 9) +
  ylim(0,60)

##POX
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(enzyme.pox, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "enzyme.pox",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "POX activity June 2021", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(enzyme.pox)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(enzyme.pox)
min(normality$p)

ggqqplot(June2021, "enzyme.pox", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
JunePOX.aov <- aov(enzyme.pox~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JunePOX.aov)
##No significance from within effects

res.aov <- anova_test(
  data = June2021, dv = enzyme.pox, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##* significance for date

p <- June2021 %>%
  pairwise_t_test(
    enzyme.pox ~ date, 
    p.adjust.method = "bonferroni"
  )

POX_June_plot<- ggline(June2021, x="date", y = "enzyme.pox", color = "treatment", title = "",
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
  geom_text(x=1, y =30000, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=30000, label = "Post", color = 'black', size = 9)+
  geom_bracket(
    xmin = c("0", "0", "0", "0"), xmax = c("1", "2", "3", "7"), 
    y.position = c(27000, 25000, 23000, 21000),
    label = c("*", "**", "*", "***"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,30000)


##Todd's Meadow
##LAP
summary <- TM2021 %>%
  group_by(sampling.date) %>%
  get_summary_stats(enzyme.lap, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "enzyme.lap",
          color = "treatment", palette = 'blue', shape = "treatment", add = c("mean_se", "jitter"),
          title = "LAP activity June 2021 Natural Pulse", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(sampling.date) %>%
  identify_outliers(enzyme.lap)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(sampling.date) %>%
  shapiro_test(enzyme.lap)
min(normality$p)

ggqqplot(TM2021, "enzyme.lap", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
TMLAP.aov <- aov(enzyme.lap~date + Error(core.id/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TMLAP.aov)
##No significance from within effects

res.aov <- anova_test(
  data = TM2021, dv = enzyme.lap, wid = core.id,
  between = date
)
get_anova_table(res.aov)
##No significance 


LAP_JuneTM_plot<- ggline(TM2021, x="date", y = "enzyme.lap", color = "treatment", title = "Todd's Meadow",
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
  geom_text(x=0.85, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=60, label = "Post", color = 'black', size = 9) + ylim (0,60) 


##AP
summary <- TM2021 %>%
  group_by(sampling.date) %>%
  get_summary_stats(enzyme.ap, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "enzyme.ap",
          color = "treatment", palette = 'blue', shape = "treatment", add = c("mean_se", "jitter"),
          title = "AP activity June Natural Pulse 2021", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(sampling.date) %>%
  identify_outliers(enzyme.ap)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(sampling.date) %>%
  shapiro_test(enzyme.ap)
min(normality$p)

ggqqplot(TM2021, "enzyme.ap", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
TMAP.aov <- aov(enzyme.ap~date + Error(core.id/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TMAP.aov)
##No significance from within effects

res.aov <- anova_test(
  data = TM2021, dv = enzyme.ap, wid = core.id,
  between = c(date)
)
get_anova_table(res.aov)
##no significance 


AP_JuneTM_plot<- ggline(TM2021, x="date", y = "enzyme.ap", color = "treatment", title = "",
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
  geom_text(x=0.85, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=60, label = "Post", color = 'black', size = 9) +
  ylim(0,60)

##June TM BG
summary <- TM2021 %>%
  group_by(sampling.date) %>%
  get_summary_stats(enzyme.bg, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "enzyme.bg",
          color = "treatment", palette = 'blue', shape = "treatment", add = c("mean_se", "jitter"),
          title = "BG activity June 2021 Natural Pulse", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(sampling.date) %>%
  identify_outliers(enzyme.bg)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(sampling.date) %>%
  shapiro_test(enzyme.bg)
min(normality$p)

ggqqplot(TM2021, "enzyme.bg", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
TMBG.aov <- aov(enzyme.bg~date + Error(core.id/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TMBG.aov)
##No significance on date

res.aov <- anova_test(
  data = TM2021, dv = enzyme.bg, wid = plot.number,
  between = c(date)
)
get_anova_table(res.aov)
##* significance

##post hoc on date effect
TM2021 %>%
  pairwise_t_test(
    enzyme.bg ~ date, 
    p.adjust.method = "bonferroni"
  )
## no sig adj.p.sig 


BG_JuneTM_plot<- ggline(TM2021, x="date", y = "enzyme.bg", color = "treatment", title = "Todd's Meadow",
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
  geom_text(x=0.85, y =60, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=60, label = "Post", color = 'black', size = 9) +
  ylim(0,60)

##POX
summary <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(enzyme.pox, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "enzyme.pox",
          color = "treatment", palette = c('blue'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "POX activity June 2021", ylab = 'Activty', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(enzyme.pox)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(enzyme.pox)
min(normality$p)

ggqqplot(TM2021, "enzyme.pox", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")

###ANOVA
TMPOX.aov <- aov(enzyme.pox~date + Error(plot.number/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TMPOX.aov)
##No significance from within effects

res.aov <- anova_test(
  data = TM2021, dv = enzyme.pox, wid = plot.number,
  between = c(date)
)
get_anova_table(res.aov)
##* significance for date

p <- TM2021 %>%
  pairwise_t_test(
    enzyme.pox ~ date, 
    p.adjust.method = "bonferroni"
  )

POX_JuneTM_plot<- ggline(TM2021, x="date", y = "enzyme.pox", color = "treatment", title = "",
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
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  labs(x = 'Day of Pulse', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =30000, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=30000, label = "Post", color = 'black', size = 9) +
  ylim(0,30000)

##Compare Sampling Campaigns

##LAP
summary <- all_pulse %>%
  group_by(sampling.campaign) %>%
  get_summary_stats(enzyme.lap, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "enzyme.lap",
          color = "treatment", shape = "treatment", add = c("mean_se", "jitter"),
          palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue'),
          title = "LAP Activity Between Sampling Campaigns", ylab = expression(LAP~(nmol~g^-1~h^-1)), 
          xlab = "Sampling Campaign", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  identify_outliers(enzyme.lap)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  shapiro_test(enzyme.lap)
min(normality$p)

ggqqplot(all_pulse, "enzyme.lap", ggtheme = theme_bw()) +
  facet_grid(. ~ sampling.campaign, labeller = "label_both")

##ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = enzyme.lap, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)
##date significance and . on sampling campaign; only want to compare sampling dates

all_pulse %>%
  pairwise_t_test(
    enzyme.lap ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )

##plot results
LAP_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "enzyme.lap",
                                   color = "black", add = c("mean_se", "jitter"),
                                   fill = 'sampling.campaign',
                                   palette =c('grey', 'orange', 'blue'),
                                   legend.title = "Sampling Campaign", legend = "bottom",
                                   title = "LAP Activity Between Sampling Campaigns", ylab = expression(LAP~(nmol~g^-1~h^-1)), 
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
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse", "GIRF September 2020 Pulse", "June 2021 GIRF Pulse"), xmax = c("June 2021 GIRF Pulse", "June 2021 TM Natural Pulse", "June 2021 TM Natural Pulse"), 
    y.position = c(66, 63, 60),
    label = c("*", "***", "***"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  ylim(0,67)


##AP
summary <- all_pulse %>%
  group_by(sampling.campaign) %>%
  get_summary_stats(enzyme.ap, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "enzyme.ap",
          color = "treatment", shape = "treatment", add = c("mean_se", "jitter"),
          palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue'),
          title = "AP Activity Between Sampling Campaigns", ylab = expression(AP~(nmol~g^-1~h^-1)), 
          xlab = "Sampling Campaign", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  identify_outliers(enzyme.ap)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  shapiro_test(enzyme.ap)
min(normality$p)

ggqqplot(all_pulse, "enzyme.ap", ggtheme = theme_bw()) +
  facet_grid(. ~ sampling.campaign, labeller = "label_both")

##ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = enzyme.ap, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)
##date significance and . on sampling campaign; only want to compare sampling dates


all_pulse %>%
  pairwise_t_test(
    enzyme.ap ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )

##plot results
AP_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "enzyme.ap",
                                  color = "black", add = c("mean_se", "jitter"),
                                  fill = 'sampling.campaign',
                                  palette =c('grey', 'orange', 'blue'),
                                  legend.title = "Sampling Campaign", legend = "bottom",
                                  title = "AP Activity Between Sampling Campaigns", ylab = expression(AP~(nmol~g^-1~h^-1)), 
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
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse", "GIRF September 2020 Pulse", "June 2021 GIRF Pulse"), xmax = c("June 2021 GIRF Pulse", "June 2021 TM Natural Pulse", "June 2021 TM Natural Pulse"), 
    y.position = c(66, 63, 60),
    label = c("**", "***", "***"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  ylim(0,67)


##BG
summary <- all_pulse %>%
  group_by( sampling.campaign) %>%
  get_summary_stats(enzyme.bg, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "enzyme.bg",
          color = "treatment", shape = "treatment", add = c("mean_se", "jitter"),
          palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue'),
          title = "BG Activity Between Sampling Campaigns", ylab = expression(BG~(nmol~g^-1~h^-1)), 
          xlab = "Sampling Campaign", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  identify_outliers(enzyme.bg)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  shapiro_test(enzyme.bg)
min(normality$p)

ggqqplot(all_pulse, "enzyme.bg", ggtheme = theme_bw()) +
  facet_grid(. ~ sampling.campaign, labeller = "label_both")


##ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = enzyme.bg, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)
##* on sampling campaign


all_pulse %>%
  pairwise_t_test(
    enzyme.bg ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )

##plot results
BG_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "enzyme.bg",
                                  color = "black", add = c("mean_se", "jitter"),
                                  fill = 'sampling.campaign',
                                  palette =c('grey', 'orange', 'blue'),
                                  legend.title = "Sampling Campaign", legend = "bottom",
                                  title = "BG Activity Between Sampling Campaigns", ylab = expression(BG~(nmol~g^-1~h^-1)), 
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
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse", "GIRF September 2020 Pulse", "June 2021 GIRF Pulse"), xmax = c("June 2021 GIRF Pulse", "June 2021 TM Natural Pulse", "June 2021 TM Natural Pulse"), 
    y.position = c(66, 63, 60),
    label = c("***", "***", "***"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  ylim(0,67)

##POX
summary <- all_pulse %>%
  group_by( sampling.campaign) %>%
  get_summary_stats(enzyme.pox, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "enzyme.pox",
          color = "treatment", shape = "treatment", add = c("mean_se", "jitter"),
          palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue'),
          title = "POX Activity Between Sampling Campaigns", ylab = expression(POX~(nmol~g^-1~h^-1)), 
          xlab = "Sampling Campaign", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  identify_outliers(enzyme.pox)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(treatment, date, sampling.campaign) %>%
  shapiro_test(enzyme.pox)
min(normality$p)

ggqqplot(all_pulse, "enzyme.pox", ggtheme = theme_bw()) +
  facet_grid(. ~ sampling.campaign, labeller = "label_both")

##ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = enzyme.pox, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)
##date significance and . on sampling campaign; only want to compare sampling dates

all_pulse %>%
  pairwise_t_test(
    enzyme.pox ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )

##plot results
POX_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "enzyme.pox",
                                   color = "black", add = c("mean_se", "jitter"),
                                   fill = 'sampling.campaign',
                                   palette =c('grey', 'orange', 'blue'),
                                   legend.title = "Sampling Campaign", legend = "bottom",
                                   title = "POX Activity Between Sampling Campaigns", ylab = expression(POX~(nmol~g^-1~h^-1)), 
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
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse", "June 2021 GIRF Pulse"), xmax = c("June 2021 GIRF Pulse","June 2021 TM Natural Pulse"), 
    y.position = c(25500, 24500),
    label = c("***", "***"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  ylim(0,26000)



# SMA for EEA -------------------------------------------------------------
##Sept 2020

##BG vs. LAP
NC_sma <- sma(lnBG~lnLAP*treatment, data=Sept2020) ##samtr major axis test if there is a common slope
NC_sma_E <- sma(lnBG~lnLAP+treatment, data=Sept2020) ##samtr major axis test for elevation (intercept) from common slope
##Accept H0: slopes are equal...p-value 0.84595, likelihood 0.3346
##Accept Ho taht there is no difference in elevation (intercepts) are not different...p-value 0.40412
NC_sma_df <- data.frame(treatment= c('Diverse', 'Grass', 'Mulch'), coef(NC_sma))


BG_LAP_2020 <- ggplot(Sept2020, aes(x=lnLAP, y=lnBG, fill = treatment)) + 
  geom_point(aes(shape = treatment, fill = treatment), size = 5)  +
  scale_fill_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  scale_shape_manual(values = c(22,24,25))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'horizontal', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  guides(color = guide_legend(order =1), shape = guide_legend(order = 2), 
         linetype= guide_legend(order = 3)) + 
  labs(x = 'lnLAP', y = 'lnBG', 
       title = 'Sept 2020 N Aquisition to C Aquisition') +
  geom_abline(aes(intercept = 0, slope = 1),size = 0.5, color = 'black', linetype = '12345678') +
  geom_abline(data=NC_sma_df,aes(intercept=elevation,slope=slope, linetype = treatment, color = treatment), size = 1.2)+
  guides(shape = guide_legend("Plot Treatment", order = 1), 
         fill = guide_legend("Plot Treatment", order = 1), 
         linetype= guide_legend(order = 2, "Line Type"),
         color = "none") + 
  labs(shape = 'Plot Treatment', fill = "Plot Treatment", linetype = "Line Type") +
  scale_linetype_manual(values = c('1:1 Reference Line' = '12345678', 'Diverse' = 'solid', "Grass" = 'dotdash', "Mulch" = "dashed")) +
  xlim(0,4) + ylim (0,4)


##BG vs. AP
NC_sma1 <- sma(lnBG~lnAP*treatment, data=Sept2020) ##samtr major axis test for common slope
NC_sma1_E <- sma(lnBG~lnAP+treatment, data=Sept2020) ##samtr major axis test for elevation (intercept) from common slope
##slopes are equal p-value 0.49095
##no difference in intercepts p-value 0.82291

NC_sma_df1 <- data.frame(treatment= c('Diverse', 'Grass', 'Mulch'), coef(NC_sma1))

BG_AP_2020 <- ggplot(Sept2020, aes(x=lnAP, y=lnBG, fill = treatment)) + 
  geom_point(aes(shape = treatment, fill = treatment), size = 5)  +
  scale_fill_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  scale_shape_manual(values = c(22,24,25))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'horizontal', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  guides(color = guide_legend(order =1), shape = guide_legend(order = 2), 
         linetype= guide_legend(order = 3)) + 
  labs(x = 'lnAP', y = 'lnBG', 
       title = 'Sept 2020 P Aquisition to C Aquisition') +
  geom_abline(aes(intercept = 0, slope = 1),size = 0.5, color = 'black', linetype = '12345678') +
  geom_abline(data=NC_sma_df1,aes(intercept=elevation,slope=slope, linetype = treatment, color = treatment), size = 1.2)+
  guides(shape = guide_legend("Plot Treatment", order = 1), 
         fill = guide_legend("Plot Treatment", order = 1), 
         linetype= guide_legend(order = 2, "Line Type"),
         color = "none") + 
  labs(shape = 'Plot Treatment', fill = "Plot Treatment", linetype = "Line Type") +
  scale_linetype_manual(values = c('1:1 Reference Line' = '12345678', 'Diverse' = 'solid', "Grass" = 'dotdash', "Mulch" = "dashed")) +
  xlim(0,4) + ylim (0,4)

##LAP vs. AP
NC_sma2 <- sma(lnLAP~lnAP*treatment, data=Sept2020) ##samtr major axis test for common slope
NC_sma2_E <- sma(lnLAP~lnAP+treatment, data=Sept2020) ##samtr major axis test for common elevation assuming equal slope
##slopes are equal p-value 0.46629
##no difference in intercept p value 0.72015

NC_sma_df2 <- data.frame(treatment= c('Diverse', 'Grass', 'Mulch'), coef(NC_sma2))

LAP_AP_2020 <- ggplot(Sept2020, aes(x=lnAP, y=lnLAP, fill = treatment)) + 
  geom_point(aes(shape = treatment, fill = treatment), size = 5)  +
  scale_fill_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  scale_shape_manual(values = c(22,24,25))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'horizontal', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  guides(color = guide_legend(order =1), shape = guide_legend(order = 2), 
         linetype= guide_legend(order = 3)) + 
  labs(x = 'lnAP', y = 'lnLAP', 
       title = 'Sept 2020 N Aquisition to P Aquisition') +
  geom_abline(aes(intercept = 0, slope = 1),size = 0.5, color = 'black', linetype = '12345678') +
  geom_abline(data=NC_sma_df2,aes(intercept=elevation,slope=slope, linetype = treatment, color = treatment), size = 1.2)+
  guides(shape = guide_legend("Plot Treatment", order = 1), 
         fill = guide_legend("Plot Treatment", order = 1), 
         linetype= guide_legend(order = 2, "Line Type"),
         color = "none") + 
  labs(shape = 'Plot Treatment', fill = "Plot Treatment", linetype = "Line Type") +
  scale_linetype_manual(values = c('1:1 Reference Line' = '12345678', 'Diverse' = 'solid', "Grass" = 'dotdash', "Mulch" = "dashed")) +
  xlim(0,4) + ylim (0,4)



##June 2021 and Todd's Meadow

##BG vs. LAP
NC_sma3 <- sma(lnBG~lnLAP*treatment, data=June2021_1) ##samtr major axis test for common slope
NC_sma3_E <- sma(lnBG~lnLAP+treatment, data=June2021_1) ##samtr major axis test for common elevation assuming equal slope
##slopes are equal p-value 0.43453
##intercepts are different p-value 0.017821


NC_sma_df3 <- data.frame(treatment= c('Diverse', 'Grass', 'Mulch', 'Reference'), coef(NC_sma3))

BG_LAP_2021 <- ggplot(June2021_1, aes(x=lnLAP, y=lnBG, fill = treatment)) + 
  geom_point(aes(shape = treatment, fill = treatment), size = 5)  +
  scale_fill_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue')) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue')) +
  scale_shape_manual(values = c(22,24,25,23))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'horizontal', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  guides(color = guide_legend(order =1), shape = guide_legend(order = 2), 
         linetype= guide_legend(order = 3)) + 
  labs(x = 'lnLAP', y = 'lnBG', 
       title = 'June 2021 N Aquisition to C Aquisition') +
  geom_abline(aes(intercept = 0, slope = 1),size = 0.5, color = 'black', linetype = '12345678') +
  geom_abline(data=NC_sma_df3,aes(intercept=elevation,slope=slope, linetype = treatment, color = treatment), size = 1.2)+
  guides(shape = guide_legend("Plot Treatment", order = 1), 
         fill = guide_legend("Plot Treatment", order = 1), 
         linetype= guide_legend(order = 2, "Line Type"),
         color = "none") + 
  labs(shape = 'Plot Treatment', fill = "Plot Treatment", linetype = "Line Type") +
  scale_linetype_manual(values = c('1:1 Reference Line' = '12345678', 'Diverse' = 'solid', "Grass" = 'dotdash', "Mulch" = "dashed", "Reference" = "dotted")) +
  xlim(0,6) + ylim (0,6)

##BG vs. AP
NC_sma4 <- sma(lnBG~lnAP*treatment, data=June2021_1) ##samtr major axis test for common slope
NC_sma4_E <- sma(lnBG~lnAP+treatment, data=June2021_1) ##samtr major axis test for common elevation assuming equal slope
##slopes are equal p-value 0.53685
##intercepts are different p-value 0.03859

NC_sma_df4 <- data.frame(treatment= c('Diverse', 'Grass', 'Mulch', 'Reference'), coef(NC_sma4))

BG_AP_2021 <- ggplot(June2021_1, aes(x=lnAP, y=lnBG, fill = treatment)) + 
  geom_point(aes(shape = treatment, fill = treatment), size = 5)  +
  scale_fill_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue')) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue')) +
  scale_shape_manual(values = c(22,24,25, 23))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'horizontal', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  guides(color = guide_legend(order =1), shape = guide_legend(order = 2), 
         linetype= guide_legend(order = 3)) + 
  labs(x = 'lnAP', y = 'lnBG', 
       title = 'June 2021 P Aquisition to C Aquisition') +
  geom_abline(aes(intercept = 0, slope = 1),size = 0.5, color = 'black', linetype = '12345678') +
  geom_abline(data=NC_sma_df4,aes(intercept=elevation,slope=slope, linetype = treatment, color = treatment), size = 1.2)+
  guides(shape = guide_legend("Plot Treatment", order = 1), 
         fill = guide_legend("Plot Treatment", order = 1), 
         linetype= guide_legend(order = 2, "Line Type"),
         color = "none") + 
  labs(shape = 'Plot Treatment', fill = "Plot Treatment", linetype = "Line Type") +
  scale_linetype_manual(values = c('1:1 Reference Line' = '12345678', 'Diverse' = 'solid', "Grass" = 'dotdash', "Mulch" = "dashed", "Reference" = "dotted")) +
  xlim(0,6) + ylim (0,6)

##LAP vs. AP
NC_sma5 <- sma(lnLAP~lnAP*treatment, data=June2021_1) ##samtr major axis test for common slope
NC_sma5_E <- sma(lnLAP~lnAP+treatment, data=June2021_1) ##samtr major axis test for common elevation assuming equal slope
##Accept Ho that the slpes are equal p-value 0.38068
##Accept Ho that there is no difference in intercept p-value 0.53217

NC_sma_df5 <- data.frame(treatment= c('Diverse', 'Grass', 'Mulch', 'Reference'), coef(NC_sma5))

LAP_AP_2021 <- ggplot(June2021_1, aes(x=lnAP, y=lnLAP, fill = treatment)) + 
  geom_point(aes(shape = treatment, fill = treatment), size = 5)  +
  scale_fill_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue')) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue')) +
  scale_shape_manual(values = c(22,24,25,23))+
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
        plot.background = element_rect((fill = 'white')),
        axis.text.x=element_text(size = 28, color = 'black'),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title.x = element_text(size = 32, color = 'black'),
        axis.title.y = element_text(size = 32, color = 'black'),
        axis.ticks=element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'))+
  theme(legend.position = 'bottom', legend.box = 'horizontal', 
        legend.background = element_rect(size = 0.5), 
        legend.title = element_text(size = 16, face = 'bold', color = 'black'))+
  guides(color = guide_legend(order =1), shape = guide_legend(order = 2), 
         linetype= guide_legend(order = 3)) + 
  labs(x = 'lnAP', y = 'lnLAP', 
       title = 'June 2021 N Aquisition to P Aquisition') +
  geom_abline(aes(intercept = 0, slope = 1),size = 0.5, color = 'black', linetype = '12345678') +
  geom_abline(data=NC_sma_df5,aes(intercept=elevation,slope=slope, linetype = treatment, color = treatment), size = 1.2)+
  guides(shape = guide_legend("Plot Treatment", order = 1), 
         fill = guide_legend("Plot Treatment", order = 1), 
         linetype= guide_legend(order = 2, "Line Type"),
         color = "none") + 
  labs(shape = 'Plot Treatment', fill = "Plot Treatment", linetype = "Line Type") +
  scale_linetype_manual(values = c('1:1 Reference Line' = '12345678', 'Diverse' = 'solid', "Grass" = 'dotdash', "Mulch" = "dashed", "Reference" = "dotted")) +
  xlim(0,6) + ylim (0,6)


##Compare Experimental Sampling Campaigns EEA SMA test
all_EEA <- filter(all_pulse, sampling.campaign != "June 2021 TM Natural Pulse")

NC_sma6 <- sma(lnLAP~lnAP*sampling.campaign, data=all_EEA) ##samtr major axis test for common slope
NC_sma6_E <- sma(lnLAP~lnAP+sampling.campaign, data=all_EEA)
##Reject Ho; slopes are not equal p-value 0.042
##Accept Ho that there is no difference in intercepts p-value 0.675

NC_sma7 <- sma(lnBG~lnAP*sampling.campaign, data=all_EEA) ##samtr major axis test for common slope
NC_sma7_E <- sma(lnBG~lnAP+sampling.campaign, data=all_EEA)
##Accept Ho that slopes are equal p-value 0.082
##Reject Ho intercepts are different p-value <0.001

NC_sma8 <- sma(lnBG~lnLAP*sampling.campaign, data=all_EEA) ##samtr major axis test for common slope
NC_sma8_E <- sma(lnBG~lnLAP+sampling.campaign, data=all_EEA)
##Accept Ho that slopes are equal p-value 0.836
##Reject Ho intercepts are different p-value <0.001



##compare Experimental:Natural Sampling Campaigns
NC_sma9 <- sma(lnLAP~lnAP*sampling.campaign, data=all_pulse) ##samtr major axis test for common slope
NC_sma9_E <- sma(lnLAP~lnAP+sampling.campaign, data=all_pulse)
##Reject Ho; slopes are not equal p-value <0.001
##Reject Ho intercepts are different p-value <0.001

NC_sma_df9 <- data.frame(sampling.campaign= c('GIRF September 2020 Pulse', 'June 2021 TM Natural Pulse', 'June 2021 GIRF Pulse'), coef(NC_sma9))

LAP_AP_sampling.campaign <- ggplot(all_pulse, aes(x=lnAP, y=lnLAP, fill = sampling.campaign)) + 
  geom_point(aes(shape = sampling.campaign, fill = sampling.campaign), size = 5)  +
  scale_fill_manual(values = c('grey', 'orange', 'blue')) +
  scale_color_manual(values = c('grey', 'orange', 'blue')) +
  scale_shape_manual(values = c(22,24,25))+
  theme_linedraw() +
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
  guides(color = guide_legend(order =1), shape = guide_legend(order = 2), 
         linetype= guide_legend(order = 3)) + 
  labs(x = 'lnAP', y = 'lnLAP', 
       title = 'Experimental:Natural Pulse N Aquisition to P Aquisition') +
  geom_abline(aes(intercept = 0, slope = 1),size = 0.5, color = 'black', linetype = '12345678') +
  geom_abline(data=NC_sma_df9,aes(intercept=elevation,slope=slope, linetype = sampling.campaign, color = sampling.campaign), size = 1.2)+
  guides(shape = guide_legend("Sampling Campaign", order = 1), 
         fill = guide_legend("Sampling Campaign", order = 1), 
         linetype= guide_legend(order = 2, "Line Type"),
         color = "none") + 
  labs(shape = 'Sampling Campaign', fill = "Sampling Campaign", linetype = "Line Type") +
  scale_linetype_manual(values = c('1:1 Reference Line' = '12345678', 'GIRF September 2020 Pulse' = 'solid', "June 2021 GIRF Pulse" = 'dotdash', "June 2021 TM Natural Pulse" = "dotted")) +
  xlim(0,6) + ylim (0,6)

NC_sma10 <- sma(lnBG~lnAP*sampling.campaign, data=all_pulse) ##samtr major axis test for common slope
NC_sma10_E <- sma(lnBG~lnAP+sampling.campaign, data=all_pulse)
##Reject Ho; slopes are not equal p-value <0.001
##Reject Ho intercepts are different p-value <0.001

NC_sma_df10 <- data.frame(sampling.campaign= c('GIRF September 2020 Pulse', 'June 2021 TM Natural Pulse', 'June 2021 GIRF Pulse'), coef(NC_sma10))

BG_AP_sampling.campaign <- ggplot(all_pulse, aes(x=lnAP, y=lnBG, fill = sampling.campaign)) + 
  geom_point(aes(shape = sampling.campaign, fill = sampling.campaign), size = 5)  +
  scale_fill_manual(values = c('grey', 'orange', 'blue')) +
  scale_color_manual(values = c('grey', 'orange', 'blue')) +
  scale_shape_manual(values = c(22,24,25))+
  theme_linedraw() +
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
  guides(color = guide_legend(order =1), shape = guide_legend(order = 2), 
         linetype= guide_legend(order = 3)) + 
  labs(x = 'lnAP', y = 'lnBG', 
       title = 'Experimental:Natural Pulse P Aquisition to C Aquisition') +
  geom_abline(aes(intercept = 0, slope = 1),size = 0.5, color = 'black', linetype = '12345678') +
  geom_abline(data=NC_sma_df10,aes(intercept=elevation,slope=slope, linetype = sampling.campaign, color = sampling.campaign), size = 1.2)+
  guides(shape = guide_legend("Sampling Campaign", order = 1), 
         fill = guide_legend("Sampling Campaign", order = 1), 
         linetype= guide_legend(order = 2, "Line Type"),
         color = "none") + 
  labs(shape = 'Sampling Campaign', fill = "Sampling Campaign", linetype = "Line Type") +
  scale_linetype_manual(values = c('1:1 Reference Line' = '12345678', 'GIRF September 2020 Pulse' = 'solid', "June 2021 GIRF Pulse" = 'dotdash', "June 2021 TM Natural Pulse" = "dotted")) +
  xlim(0,6) + ylim (0,6)

NC_sma11 <- sma(lnBG~lnLAP*sampling.campaign, data=all_pulse) ##samtr major axis test for common slope
NC_sma11_E <- sma(lnBG~lnLAP+sampling.campaign, data=all_pulse)
##Accept Ho that slopes are equal p-value 0.150
##Reject Ho intercepts are different p-value <0.001

NC_sma_df11 <- data.frame(sampling.campaign= c('GIRF September 2020 Pulse', 'June 2021 TM Natural Pulse', 'June 2021 GIRF Pulse'), coef(NC_sma11))

BG_LAP_sampling.campaign <- ggplot(all_pulse, aes(x=lnLAP, y=lnBG, fill = sampling.campaign)) + 
  geom_point(aes(shape = sampling.campaign, fill = sampling.campaign), size = 5)  +
  scale_fill_manual(values = c('grey', 'orange', 'blue')) +
  scale_color_manual(values = c('grey', 'orange', 'blue')) +
  scale_shape_manual(values = c(22,24,25))+
  theme_linedraw() +
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
  guides(color = guide_legend(order =1), shape = guide_legend(order = 2), 
         linetype= guide_legend(order = 3)) + 
  labs(x = 'lnLAP', y = 'lnBG', 
       title = 'Experimental:Natural Pulse N Aquisition to C Aquisition') +
  geom_abline(aes(intercept = 0, slope = 1),size = 0.5, color = 'black', linetype = '12345678') +
  geom_abline(data=NC_sma_df11,aes(intercept=elevation,slope=slope, linetype = sampling.campaign, color = sampling.campaign), size = 1.2)+
  guides(shape = guide_legend("Sampling Campaign", order = 1), 
         fill = guide_legend("Sampling Campaign", order = 1), 
         linetype= guide_legend(order = 2, "Line Type"),
         color = "none") + 
  labs(shape = 'Sampling Campaign', fill = "Sampling Campaign", linetype = "Line Type") +
  scale_linetype_manual(values = c('1:1 Reference Line' = '12345678', 'GIRF September 2020 Pulse' = 'solid', "June 2021 GIRF Pulse" = 'dotdash', "June 2021 TM Natural Pulse" = "dotted")) +
  xlim(0,6) + ylim (0,6)


##SMA post hoc
JuneG_D <- filter(June2021, treatment == "Grass" | treatment == "Diverse")
JuneG_M <- filter(June2021, treatment == "Grass" | treatment == "Mulch")
JuneD_M <- filter(June2021, treatment == "Diverse" | treatment == "Mulch")


inter_JuneGD <- sma(lnBG~lnLAP+treatment, data=JuneG_D)
inter_JuneGM <- sma(lnBG~lnLAP+treatment, data=JuneG_M)
inter_JuneDM <- sma(lnBG~lnLAP+treatment, data=JuneD_M)

inter_JuneGD1 <- sma(lnBG~lnAP+treatment, data=JuneG_D)
inter_JuneGM1 <- sma(lnBG~lnAP+treatment, data=JuneG_M)
inter_JuneDM1 <- sma(lnBG~lnAP+treatment, data=JuneD_M)

JT <- filter(all_pulse, sampling.campaign == "June 2021 GIRF Pulse" | sampling.campaign == "June 2021 TM Natural Pulse")
ST <- filter(all_pulse, sampling.campaign == "GIRF September 2020 Pulse" | sampling.campaign == "June 2021 TM Natural Pulse")
SJ <- filter(all_pulse, sampling.campaign == "GIRF September 2020 Pulse" | sampling.campaign == "June 2021 GIRF Pulse")

slope_JT <- sma(lnBG~lnAP*sampling.campaign, data=JT)
slope_ST <- sma(lnBG~lnAP*sampling.campaign, data=ST)
slope_SJ <- sma(lnBG~lnAP*sampling.campaign, data=SJ)

slope_JT1 <- sma(lnLAP~lnAP*sampling.campaign, data=JT)
slope_ST1 <- sma(lnLAP~lnAP*sampling.campaign, data=ST)
slope_SJ1 <- sma(lnLAP~lnAP*sampling.campaign, data=SJ)

inter_JT <- sma(lnBG~lnAP+sampling.campaign, data=JT)
inter_ST <- sma(lnBG~lnAP+sampling.campaign, data=ST)
inter_SJ <- sma(lnBG~lnAP+sampling.campaign, data=SJ)

inter_JT1 <- sma(lnLAP~lnAP+sampling.campaign, data=JT)
inter_ST1 <- sma(lnLAP~lnAP+sampling.campaign, data=ST)
inter_SJ1 <- sma(lnLAP~lnAP+sampling.campaign, data=SJ)

inter_JT2 <- sma(lnBG~lnLAP+sampling.campaign, data=JT)
inter_ST2 <- sma(lnBG~lnLAP+sampling.campaign, data=ST)
inter_SJ2 <- sma(lnBG~lnLAP+sampling.campaign, data=SJ)



# Vector Analysis ---------------------------------------------------------
summary <- Sept2020 %>%
  get_summary_stats(length | angle, type = 'mean_se') 
summary

summary1 <- June2021 %>%
  get_summary_stats(length | angle, type = 'mean_se') 
summary1

summary2 <- TM2021 %>%
  get_summary_stats(length | angle, type = 'mean_se') 
summary2

##ANOVA on length
res.aov <- anova_test(
  data = all_pulse, dv = length, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)
##*** significance

pwc <- all_pulse %>%
  pairwise_t_test(
    length ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )

all_length <- ggboxplot(all_pulse, x = "sampling.campaign", y = "length",
                        fill = 'sampling.campaign',
                        palette =c('grey', 'orange', 'blue'),
                        legend.title = "Sampling Campaign", legend = "bottom", add = c("mean_se", "jitter"),
                        title = "Vector Length Between Pulse Events", ylab = 'Vector Length', 
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
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse", "GIRF September 2020 Pulse"), xmax = c("June 2021 GIRF Pulse","June 2021 TM Natural Pulse"),
    y.position = c(2.8, 2.6),
    label = c("***", "***"),
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  scale_y_continuous(limits = c(0,3))

##ANOVA on angle
res.aov <- anova_test(
  data = all_pulse, dv = angle, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)
##* significance

pwc <- all_pulse %>%
  pairwise_t_test(
    angle ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )

all_angle  <- ggboxplot(all_pulse, x = "sampling.campaign", y = "angle",
                        fill = 'sampling.campaign',
                        palette =c('grey', 'orange', 'blue'),
                        legend.title = "Sampling Campaign", legend = "bottom", add = c("mean_se", "jitter"),
                        title = "Vector Angle Between Pulse Events", ylab = 'Vector Angle (?)', 
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
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse"), xmax = c("June 2021 GIRF Pulse"),
    y.position = c(80),
    label = c("*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  geom_hline(yintercept=45, linetype="dashed", 
             color = "grey48", size=1) +
  scale_y_continuous(limits = c(0,90)) 

##Ratio plot of sampling campaigns
vector_sampling.campaign <- ggplot(all_pulse, aes(x=lnBG_lnAP, y=lnBG_lnLAP, fill = sampling.campaign)) +
  geom_point(aes(shape = sampling.campaign, fill = sampling.campaign), size = 5)  +
  scale_fill_manual(values = c('grey', 'orange', 'blue')) +
  scale_color_manual(values = c('grey', 'orange', 'blue')) +
  scale_shape_manual(values = c(22,24,25))+
  theme_classic() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = 'black', size = 1), 
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
  guides(color = guide_legend(order =1), shape = guide_legend(order = 2),
         linetype= guide_legend(order = 3)) +
  labs(x = 'lnBG/lnAP', y = 'lnBG/lnLAP',
       title = '') +
  geom_abline(aes(intercept = 0, slope = 1),size = 1, color = 'black', linetype = 'solid') +
  guides(shape = guide_legend("Sampling Campaign", order = 1),
         fill = guide_legend("Sampling Campaign", order = 1),
         linetype= guide_legend(order = 2, "Line Type"),
         color = "none") +
  labs(shape = 'Sampling Campaign', fill = "Sampling Campaign") +
  xlim(0,2) + ylim (0,2)



ggarrange(LAP_Sept_plot, LAP_June_plot, LAP_JuneTM_plot,
          AP_Sept_plot, AP_June_plot, AP_JuneTM_plot,
          ncol = 3, nrow = 2, legend = 'none', labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
          font.label = list(size = 30))

ggarrange(BG_Sept_plot, BG_June_plot, BG_JuneTM_plot,
          POX_Sept_plot, POX_June_plot, POX_JuneTM_plot,
          ncol = 3, nrow = 2, legend = 'none', labels = c('A', 'B', 'C', 'D', 'E', 'F'), 
          font.label = list(size = 30))

ggarrange(LAP_sampling.campaign, AP_sampling.campaign,
          BG_sampling.campaign, POX_sampling.campaign,
          ncol = 2, nrow = 2, legend = 'bottom', common.legend = TRUE, labels = c('A', 'B', 'C', 'D'), 
          font.label = list(size = 30))

ggarrange(BG_LAP_2020, BG_LAP_2021, BG_LAP_sampling.campaign,
          BG_AP_2020, BG_AP_2021, BG_AP_sampling.campaign,
          LAP_AP_2020, LAP_AP_2021, LAP_AP_sampling.campaign,
          ncol = 3, nrow = 3, legend = 'none', labels = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'), 
          font.label = list(size = 30))
