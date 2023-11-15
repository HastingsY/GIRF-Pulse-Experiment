##Gravimetric Soil Moisture and pH
##Yvette Hastings
##April 17, 2022
##Updated 9/14/2023

## Data analysis sections; quick access using shift + alt + J
## 1. Data clean-up
## 2. Set ggplot theme
## 3. Gravimetric Soil Moisture
## 4. pH
## 5. Final Legends
## 6. Final plots


##load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(theme)
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

##format data
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



# Set ggplot theme --------------------------------------------------------
plot_theme <- function(){
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), axis.line = element_line(colour = '#7F7F7F', size = 1), 
    plot.background = element_rect((fill = 'white')),
    axis.text.x=element_text(size = 28, color = 'black'),
    axis.text.y = element_text(size = 28, color = "black"),
    axis.title.x = element_text(size = 28, color = 'black'),
    axis.title.y = element_text(size = 32, color = 'black'),
    axis.ticks=element_blank(),
    plot.title = element_text(hjust = 0.5, face = 'bold', size = 32, color = 'black'),
    legend.position = 'bottom', legend.box = 'vertical', 
    legend.background = element_rect(size = 0.5), 
    legend.title = element_text(size = 30, face = 'bold', color = 'black'),
    legend.key.size = unit(3, 'cm'), legend.text = element_text(size = 30)
  )
}



# Gravimetric Soil Moisture -----------------------------------------------

##Sept 2020
summary <- Sept2020 %>%
  group_by(sampling.date) %>%
  get_summary_stats(gravimetric.sm, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "gravimetric.sm",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "Gravimetric Moisture by Treatment over Sept 2020 Pulse Experiment", ylab = 'Gravimetric Moisture (%)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(gravimetric.sm)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(gravimetric.sm)
min(normality$p)

ggqqplot(Sept2020, "gravimetric.sm", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")
##the points fall along the line in the qq plot

###ANOVA
SeptM.aov <- aov(gravimetric.sm~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(SeptM.aov)
##Treatment Effect

res.aov <- anova_test(
  data = Sept2020, dv = gravimetric.sm, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##Between effect significance by treatment and date, no interaction terms

##post hoc on date effect
Sept2020 %>%
  pairwise_t_test(
    gravimetric.sm ~ date, 
    p.adjust.method = "bonferroni"
  )

##post hoc on treatment effect
Sept2020 %>%
  pairwise_t_test(
    gravimetric.sm ~ treatment, 
    p.adjust.method = "bonferroni"
  )
##Mulch shows to have the biggest influence on moisutre

Sept_moisture <- ggline(Sept2020, x="date", y = "gravimetric.sm", color = "treatment", title = "September 2020",
                        add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                        size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  plot_theme() +
  labs(x = 'Day of Pulse', y= 'Gravimetric Soil Moisture (%)') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =40, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=40, label = "Post", color = 'black', size = 9) +
  geom_bracket(
    xmin = c("0"), xmax = c("3"),
    y.position = c(25),
    label = c("*"),label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,40)

Sept_moisture_treatment <- ggboxplot(Sept2020, x = "treatment", y = "gravimetric.sm",
                                     color = "black", add = c("mean_se", "jitter"),
                                     fill = 'treatment',
                                     palette =c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                                     title = "", ylab = 'Gravimetric Moisture (%)', 
                                     xlab = "", bxp.errorbar = TRUE)+
  plot_theme() +
  geom_bracket(
    xmin = c("Diverse","Grass"), xmax = c("Mulch", "Mulch"),
    y.position = c(25, 22),
    label = c("***","***"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,30)


##June 2021
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(gravimetric.sm, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "gravimetric.sm",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "Gravimetric Moisture by Treatment over June 2021 Pulse Experiment", ylab = 'Gravimetric Moisture (%)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(gravimetric.sm)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(gravimetric.sm)
min(normality$p)

ggqqplot(June2021, "gravimetric.sm", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")
#majority of points fall on line

###ANOVA
JuneM.aov <- aov(gravimetric.sm~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(JuneM.aov)
##Date Effect

res.aov <- anova_test(
  data = June2021, dv = gravimetric.sm, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##Between effect significance by treatment and date, no interaction terms

##post hoc to show date effect
p <- June2021 %>%
  pairwise_t_test(
    gravimetric.sm ~ date, 
    p.adjust.method = "bonferroni"
  )

##post hoc on treatment effect
June2021 %>%
  pairwise_t_test(
    gravimetric.sm ~ treatment, 
    p.adjust.method = "bonferroni"
  )
##Mulch shows to have the biggest influence on moisutre

##post hoc on interaction effect
pwc <- June2021 %>%
  group_by(date) %>%
  pairwise_t_test(
    gravimetric.sm ~ treatment, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc
##Mulch shows to have the biggest difference in moisutre

##June moisture plot
June_moisture <- ggline(June2021, x="date", y = "gravimetric.sm", color = "treatment", title = "June 2021",
                        add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                        size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  plot_theme() +
  labs(x = 'Day of Pulse', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =40, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=40, label = "Post", color = 'black', size = 9)+ 
  geom_bracket(
    xmin = c("0","0","0", "0", "0", "1"), xmax = c("1", "2","3", "7", "9", "16"),
    y.position = c(35,33, 31,29,27, 25),
    label = c("***","***","***", "***", "*", "*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,40)

##treatment plot
June_moisture_treatment <- ggboxplot(June2021, x = "treatment", y = "gravimetric.sm",
                                     color = "black", add = c("mean_se", "jitter"),
                                     fill = 'treatment',
                                     palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'),
                                     title = "", ylab = '', 
                                     xlab = "", bxp.errorbar = TRUE)+
  plot_theme() +
  geom_bracket(
    xmin = c("Diverse","Grass"), xmax = c("Mulch", "Mulch"),
    y.position = c(30, 28),
    label = c("***","***"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,30)


##Todd's Meadow
summary <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(gravimetric.sm, type = 'full')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "gravimetric.sm",
          color = "treatment", palette = 'blue', shape = "treatment", add = c("mean_se", "jitter"),
          title = "Gravimetric Moisture over June 2021 Natural Pulse Experiment", ylab = 'Gravimetric Moisture (%)', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(sampling.date) %>%
  identify_outliers(gravimetric.sm)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(sampling.date) %>%
  shapiro_test(gravimetric.sm)
min(normality$p)

ggqqplot(TM2021, "gravimetric.sm", ggtheme = theme_bw()) 
##majority of points fall along the line

###ANOVA
TM_M.aov <- aov(gravimetric.sm~date + Error(core.id/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TM_M.aov)
##No significance

res.aov <- anova_test(
  data = TM2021, dv = gravimetric.sm, wid = core.id,
  between = c(date)
)
get_anova_table(res.aov)
##Between effect significance by date

##post hoc to show date effect
TM2021 %>%
  pairwise_t_test(
    gravimetric.sm ~ date, 
    p.adjust.method = "bonferroni"
  )

##June moisture plot
TM_moisture <- ggline(TM2021, x="date", y = "gravimetric.sm", color = "treatment", title = "Todd's Meadow",
                      add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                      size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = 'blue') +
  plot_theme() +
  labs(x = 'Day of Pulse', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =40, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=40, label = "Post", color = 'black', size = 9) +
  geom_bracket(
    xmin = c("0","1"), xmax = c("3", "3"),
    y.position = c(31,29),
    label = c("***","***"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)+ ylim (0,40)

TM_moisture_treatment <- ggboxplot(TM2021, x = "treatment", y = "gravimetric.sm",
                                   color = "black", add = c("mean_se", "jitter"),
                                   fill = 'treatment',
                                   palette = c('blue'),
                                   title = "", ylab = '', 
                                   xlab = "", bxp.errorbar = TRUE, width = 0.26)+
  plot_theme() +
  ylim(0,30)



##Between Sampling Campaigns
summary <- all_pulse %>%
  group_by(sampling.campaign) %>%
  get_summary_stats(gravimetric.sm, type = 'mean_se')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "gravimetric.sm",
          color = "sampling.campaign", palette = c('grey', 'orange', 'blue'), shape = "sampling.campaign", add = c("mean_se", "jitter"),
          title = "GM between Experimental and Natural Pulse", ylab = 'GM', 
          xlab = "", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(sampling.campaign) %>%
  identify_outliers(gravimetric.sm)
outliers
##outliers detected but not removed because individual tests do not show outliers

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(sampling.campaign) %>%
  shapiro_test(gravimetric.sm)
min(normality$p)

ggqqplot(all_pulse, "gravimetric.sm", ggtheme = theme_bw()) +
  facet_grid(~sampling.campaign, labeller = "label_both")
##majority of points fall along the line

###ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = gravimetric.sm, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)
##sampling campaign signifcance

all_pulse %>%
  pairwise_t_test(
    gravimetric.sm ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )

GM_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "gravimetric.sm",
                                  color = "black", add = c("mean_se", "jitter"),fill = 'sampling.campaign',
                                  palette =c('grey', 'orange', 'blue'),
                                  title = "Gravimetric Soil Moisture Between Pulse Events", ylab = 'Gravimetric Soil Moisture (%)', 
                                  xlab = "", bxp.errorbar = TRUE, legend.title = "Sampling Campaign", legend = "bottom")+
  plot_theme() +
  scale_x_discrete(labels=c("GIRF September 2020 Pulse" = "GIRF Sept 2020", "June 2021 GIRF Pulse" = "GIRF June 2021", "June 2021 TM Natural Pulse" = "Reference June 2021")) +
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse", "GIRF September 2020 Pulse"), xmax = c("June 2021 GIRF Pulse", "June 2021 TM Natural Pulse"), 
    y.position = c(34, 32),
    label = c("***", "**"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  ylim(0,35)



# pH ----------------------------------------------------------------------

##Sept 2020
summary <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(ph, type = 'mean_se')
summary

##check for outliers
ggboxplot(Sept2020, x = "date", y = "ph",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "pH by Treatment over Sept 2020 Pulse Experiment", ylab = 'pH', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(ph)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- Sept2020 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(ph)
min(normality$p)
##p>0.05, sample size is >30, continue with anova 

ggqqplot(Sept2020, "ph", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")
##the points fall along the line in the qq plot

###ANOVA
Septph.aov <- aov(ph~treatment*date + Error(plot.number/(treatment*date)), data = Sept2020) ##check interaction terms for sampling date and treatment
summary(Septph.aov)
##No significance

res.aov <- anova_test(
  data = Sept2020, dv = ph, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##Significance by date

##post hoc on date effect
Sept2020 %>%
  pairwise_t_test(
    ph ~ date, 
    p.adjust.method = "bonferroni"
  )


Sept_ph <- ggline(Sept2020, x="date", y = "ph", color = "treatment", title = "September 2020",
                  add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                  size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  plot_theme() +
  labs(x = 'Day of Pulse', y= 'pH') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =10, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=10, label = "Post", color = 'black', size = 9) +
  ylim(7,10) +
  geom_bracket(
    xmin = c("0", "3", "5", "5"), xmax = c("5", "9", "9", "16"),
    y.position = c(9.75,9.5,9.2,8.9),
    label = c("***", "***", "***", "*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)

##June 2021
summary <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(ph, type = 'mean_se')
summary

##check for outliers
ggboxplot(June2021, x = "date", y = "ph",
          color = "treatment", palette = c('darkgoldenrod1', 'forestgreen', 'chocolate4'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "pH by Treatment over June 2021 Pulse Experiment", ylab = 'pH', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(ph)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- June2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(ph)
min(normality$p)

ggqqplot(June2021, "ph", ggtheme = theme_bw()) +
  facet_grid(. ~ treatment, labeller = "label_both")
##the points fall along the line in the qq plot

###ANOVA
Juneph.aov <- aov(ph~treatment*date + Error(plot.number/(treatment*date)), data = June2021) ##check interaction terms for sampling date and treatment
summary(Juneph.aov)
##No significance

res.aov <- anova_test(
  data = June2021, dv = ph, wid = plot.number,
  between = c(treatment, date)
)
get_anova_table(res.aov)
##Significance by date

##post hoc on date effect
June2021 %>%
  pairwise_t_test(
    ph ~ date, 
    p.adjust.method = "bonferroni"
  )


June_ph <- ggline(June2021, x="date", y = "ph", color = "treatment", title = "June 2021",
                  add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                  size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4')) +
  plot_theme()+
  labs(x = 'Day of Pulse', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.6, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=1, y =10, label = "Pre", color = 'black', size = 9) +
  geom_text(x=2.3, y=10, label = "Post", color = 'black', size = 9)+ 
  ylim(7,10) +
  geom_bracket(
    xmin = c("1"), xmax = c("9"),
    y.position = c(9.5),
    label = c("*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)

##Todd's Meadow
summary <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  get_summary_stats(ph, type = 'mean_se')
summary

##check for outliers
ggboxplot(TM2021, x = "date", y = "ph",
          color = "treatment", palette = c('blue'), shape = "treatment", add = c("mean_se", "jitter"),
          title = "pH by Treatment over TM 2021 Natural Pulse", ylab = 'pH', 
          xlab = "Day of Pulse", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  identify_outliers(ph)
outliers
##no extreme outliers detected

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- TM2021 %>%
  group_by(treatment, sampling.date) %>%
  shapiro_test(ph)
min(normality$p)

ggqqplot(TM2021, "ph", ggtheme = theme_bw()) 
##the points fall along the line in the qq plot

###ANOVA
TMph.aov <- aov(ph~date + Error(plot.number/(date)), data = TM2021) ##check interaction terms for sampling date and treatment
summary(TMph.aov)
##date significance

res.aov <- anova_test(
  data = TM2021, dv = ph, wid = plot.number,
  between =  date
)
get_anova_table(res.aov)
##Significance by date

##post hoc on date effect
TM2021 %>%
  pairwise_t_test(
    ph ~ date, 
    p.adjust.method = "bonferroni"
  )


TM_ph <- ggline(TM2021, x="date", y = "ph", color = "treatment", title = "Todd's Meadow",
                add = c("mean_se", "jitter"), legend.title = "Treatment", legend = 'right', 
                size =1, shape = "treatment", point.size = 5) +
  scale_color_manual(values = c('blue')) +
  plot_theme() +
  labs(x = 'Day of Pulse', y= '') +
  labs(shape = 'Plot Treatment', color = "Plot Treatment")  +
  geom_vline(xintercept = 1.3, linetype = 'dotted', color = 'black', size = 1.5) +
  geom_text(x=0.85, y =10, label = "Pre", color = 'black', size = 9) +
  geom_text(x=1.85, y=10, label = "Post", color = 'black', size = 9) +
  ylim(7,10) +
  geom_bracket(
    xmin = c("0"), xmax = c("3"),
    y.position = c(8.5),
    label = c("*"),
    label.size = 10,
    tip.length = 0.01, size = 1.5)


##Between Sampling Campaigns
summary <- all_pulse %>%
  #group_by(sampling.campaign) %>%
  get_summary_stats(ph, type = 'full')
summary

##check for outliers
ggboxplot(all_pulse, x = "sampling.campaign", y = "ph",
          color = "sampling.campaign", palette = c('grey', 'orange', 'blue'), shape = "sampling.campaign", add = c("mean_se", "jitter"),
          title = "pH between Experimental and Natural Pulse", ylab = 'pH', 
          xlab = "", bxp.errorbar = TRUE)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  theme(legend.position = 'right', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5, colour = 'black'), 
        legend.title = element_text(size = 12, face = 'bold'))

outliers <- all_pulse %>%
  group_by(sampling.campaign) %>%
  identify_outliers(ph)
outliers
##outliers detected but not removed

##normality check
##Shapiro-Wilk's test since dataset is small
normality <- all_pulse %>%
  group_by(sampling.campaign) %>%
  shapiro_test(ph)
min(normality$p)

ggqqplot(all_pulse, "ph", ggtheme = theme_bw()) +
  facet_grid(~sampling.campaign, labeller = "label_both")
##majority of points fall along the line

###ANOVA
res.aov <- anova_test(
  data = all_pulse, dv = ph, wid = plot.number,
  between = c(sampling.campaign)
)
get_anova_table(res.aov)
##sampling campaign signifcance

all_pulse %>%
  pairwise_t_test(
    ph ~ sampling.campaign, 
    p.adjust.method = "bonferroni"
  )


ph_sampling.campaign <- ggboxplot(all_pulse, x = "sampling.campaign", y = "ph",
                                  color = "black", add = c("mean_se", "jitter"), fill = 'sampling.campaign',
                                  palette =c('grey', 'orange', 'blue'),
                                  title = "pH Between Pulse Events", ylab = 'pH', 
                                  xlab = "", bxp.errorbar = TRUE, 
                                  legend.title = "Sampling Campaign", legend = 'bottom')+
  plot_theme()+
  scale_x_discrete(labels=c("GIRF September 2020 Pulse" = "GIRF Sept 2020", "June 2021 GIRF Pulse" = "GIRF June 2021", "June 2021 TM Natural Pulse" = "Reference June 2021")) +
  geom_bracket(
    xmin = c("GIRF September 2020 Pulse", "June 2021 GIRF Pulse"), xmax = c("June 2021 GIRF Pulse", "June 2021 TM Natural Pulse"), 
    y.position = c(10, 9.8),
    label = c("***", "***"), 
    label.size = 10,
    tip.length = 0.01, size = 1.5) +
  ylim(7,10)


# Final Legends -----------------------------------------------------------
##this code is just to extract and add combined legends 
legend_1 <- ggline(all_pulse, x="date", y = "ph", color = "treatment",  legend = 'bottom', 
                size =1, shape = "treatment", point.size = 8.5, legend.title = "Treatment") +
  scale_color_manual(values = c('darkgoldenrod1', 'forestgreen', 'chocolate4', 'blue')) +
  scale_shape_manual(values = c(19, 17, 15, 19)) +
  theme(legend.position = 'bottom', legend.box = 'vertical', 
        legend.background = element_rect(size = 0.5),  
        legend.title = element_text(size = 30, face = 'bold', color = 'black'),
        legend.key.size = unit(3, 'cm'), legend.text = element_text(size = 30))

legend1 <- cowplot::get_legend(legend_1)
legend2 <- cowplot::get_legend(GM_sampling.campaign)
combined_legends <- cowplot::plot_grid(legend1, legend2, align = "hv", nrow = 1,
                                       rel_widths = c(0.55, 0.85))


# Final Plots -------------------------------------------------------------
moisture <- cowplot::plot_grid(Sept_moisture + theme(legend.position="none"), 
                               June_moisture + theme(legend.position="none"), 
                               TM_moisture + theme(legend.position="none"), 
                               Sept_moisture_treatment + theme(legend.position="none"), 
                               June_moisture_treatment + theme(legend.position="none"), 
                               TM_moisture_treatment + theme(legend.position="none"),
                               ncol = 3, nrow = 2, labels = "AUTO", label_size = 30)
moisutre_legend <- cowplot::plot_grid(moisture,
                                      legend1,
                                      ncol = 1, rel_heights = c(0.95, 0.05))


ggarrange(GM_sampling.campaign, ph_sampling.campaign,
          ncol = 2, nrow = 1, common.legend = TRUE, legend = 'bottom', labels = c('A', 'B'), 
          font.label = list(size = 30))


pH <- cowplot::plot_grid(Sept_ph + theme(legend.position="none"),
                         June_ph + theme(legend.position="none"),
                         TM_ph + theme(legend.position="none"),
                         nrow = 1, ncol = 3, labels = "AUTO", label_size = 30)
pH_legend <- cowplot::plot_grid(pH,
                                legend1,
                                ncol = 1, rel_heights = c(0.9,0.1))




