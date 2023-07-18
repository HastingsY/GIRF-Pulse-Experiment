##GIRF 2020-2021 and TM 2021 WEO Soil Sensor
##Authored by: Yvette  D. Hastings
##Written: May 25, 2022
##Updated: June 2, 2023


library(readxl)
library(dplyr)
library(dtplyr)
library(sqldf)
library(lubridate)
library(rstatix)
library(tidyverse)

##2020
P1_2020 <- read.csv("Data/dat/2020/P1.csv") 
P2_2020 <- read.csv("Data/dat/2020/P2.csv")
P3_2020 <- read.csv("Data/dat/2020/P3.csv")
P4_2020 <- read.csv("Data/dat/2020/P4.csv")
P5_2020 <- read.csv("Data/dat/2020/P5.csv")
P6_2020 <- read.csv("Data/dat/2020/P6.csv")
P7_2020 <- read.csv("Data/dat/2020/P7.csv")
P8_2020 <- read.csv("Data/dat/2020/P8.csv")
P9_2020 <- read.csv("Data/dat/2020/P9.csv")


##2021
P1_2021 <- read.csv("Data/dat/2021/P1.csv")
P2_2021 <- read.csv("Data/dat/2021/P2.csv") 
P3_2021 <- read.csv("Data/dat/2021/P3.csv") 
P4_2021 <- read.csv("Data/dat/2021/P4.csv") 
P5_2021 <- read.csv("Data/dat/2021/P5.csv") 
P6_2021 <- read.csv("Data/dat/2021/P6.csv") 
P7_2021 <- read.csv("Data/dat/2021/P7.csv") 
P8_2021 <- read.csv("Data/dat/2021/P8.csv") 
P9_2021 <- read.csv("Data/dat/2021/P9.csv") 

##remove NA and filter
P1_2020 <- P1_2020 %>%
  filter(LocalDateTime > '5/1/2020') %>%
  filter(LocalDateTime < '9/30/2020') %>%
  subset(VWC_20cm > 0 & VWC_20cm < 90) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P2_2020 <- P2_2020 %>%
  filter(LocalDateTime > '5/1/2020') %>%
  filter(LocalDateTime < '9/30/2020') %>%
  subset(VWC_20cm > 0 & VWC_20cm < 90) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P3_2020 <- P3_2020 %>%
  filter(LocalDateTime > '5/1/2020') %>%
  filter(LocalDateTime < '9/30/2020') %>%
  subset(VWC_20cm > 0 & VWC_20cm < 90) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P4_2020 <- P4_2020 %>%
  filter(LocalDateTime > '5/1/2020') %>%
  filter(LocalDateTime < '9/30/2020') %>%
  subset(VWC_20cm > 0 & VWC_20cm < 90) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P5_2020 <- P5_2020 %>%
  filter(LocalDateTime > '5/1/2020') %>%
  filter(LocalDateTime < '9/30/2020') %>%
  subset(VWC_20cm > 0 & VWC_20cm < 90) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P6_2020 <- P6_2020 %>%
  filter(LocalDateTime > '5/1/2020') %>%
  filter(LocalDateTime < '9/30/2020') %>%
  subset(VWC_20cm > 0 & VWC_20cm < 90) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P7_2020 <- P7_2020 %>%
  filter(LocalDateTime > '5/1/2020') %>%
  filter(LocalDateTime < '9/30/2020') %>%
  subset(VWC_20cm > 0 & VWC_20cm < 90) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P8_2020 <- P8_2020 %>%
  filter(LocalDateTime > '5/1/2020') %>%
  filter(LocalDateTime < '9/30/2020') %>%
  subset(VWC_20cm > 0 & VWC_20cm < 90) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P9_2020 <- P9_2020 %>%
  filter(LocalDateTime > '5/1/2020') %>%
  filter(LocalDateTime < '9/30/2020') %>%
  subset(VWC_20cm > 0 & VWC_20cm < 90) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

##join GIRF 2020 SM
GIRF2020_list <- list(P1_2020, P2_2020, P3_2020, P4_2020, P5_2020, P6_2020, P7_2020, P8_2020, P9_2020)
GIRF2020_SM_sensor <- GIRF2020_list %>%
  reduce(full_join, by = 'LocalDateTime') %>%
  rowwise() %>%
  mutate(average = mean(c_across(starts_with("VWC")), na.rm=TRUE)) %>%
  ungroup() %>%
  get_summary_stats(average, type = 'full') 

GIRF2020_SM_sensor_oneweekbeforepulse <- GIRF2020_list %>%
  reduce(full_join, by = 'LocalDateTime') %>%
  mutate(LocalDateTime = as.Date(as.character(LocalDateTime), format = "%m/%d/%Y")) %>%
  filter(LocalDateTime > '2020/08/27' & LocalDateTime < '2020/09/03') %>%
  rowwise() %>%
  mutate(average = mean(c_across(starts_with("VWC")), na.rm=TRUE)) %>%
  ungroup() %>%
  get_summary_stats(average, type = 'full') 

GIRF2020_SM_sensor_pre_pulse <- GIRF2020_list %>%
  reduce(full_join, by = 'LocalDateTime') %>%
  mutate(LocalDateTime = as.Date(as.character(LocalDateTime), format = "%m/%d/%Y")) %>%
  filter(LocalDateTime == '2020/09/03') %>%
  rowwise() %>%
  mutate(average = mean(c_across(starts_with("VWC")), na.rm=TRUE)) %>%
  ungroup() %>%
  get_summary_stats(average, type = 'full') 

GIRF2020_SM_sensor_pulse <- GIRF2020_list %>%
  reduce(full_join, by = 'LocalDateTime') %>%
  mutate(LocalDateTime = as.Date(as.character(LocalDateTime), format = "%m/%d/%Y")) %>%
  filter(LocalDateTime > '2020/09/03' & LocalDateTime < '2020/09/23') %>%
  rowwise() %>%
  mutate(average = mean(c_across(starts_with("VWC")), na.rm=TRUE)) %>%
  ungroup() %>%
  get_summary_stats(average, type = 'full') 


P1_2021 <- P1_2021 %>%
  filter(LocalDateTime > '5/1/2021') %>%
  filter(LocalDateTime < '9/30/2021') %>%
  subset(VWC_20cm > 0) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P2_2021 <- P2_2021 %>%
  filter(LocalDateTime > '5/1/2021') %>%
  filter(LocalDateTime < '9/30/2021') %>%
  subset(VWC_20cm > 0) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P3_2021 <- P3_2021 %>%
  filter(LocalDateTime > '5/1/2021') %>%
  filter(LocalDateTime < '9/30/2021') %>%
  subset(VWC_20cm > 0) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P4_2021 <- P4_2021 %>%
  filter(LocalDateTime > '5/1/2021') %>%
  filter(LocalDateTime < '9/30/2021') %>%
  subset(VWC_20cm > 0) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P5_2021 <- P5_2021 %>%
  filter(LocalDateTime > '5/1/2021') %>%
  filter(LocalDateTime < '9/30/2021') %>%
  subset(VWC_20cm > 0) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P6_2021 <- P6_2021 %>%
  filter(LocalDateTime > '5/1/2021') %>%
  filter(LocalDateTime < '9/30/2021') %>%
  subset(VWC_20cm > 0) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P7_2021 <- P7_2021 %>%
  filter(LocalDateTime > '5/1/2021') %>%
  filter(LocalDateTime < '9/30/2021') %>%
  subset(VWC_20cm > 0) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P8_2021 <- P8_2021 %>%
  filter(LocalDateTime > '5/1/2021') %>%
  filter(LocalDateTime < '9/30/2021') %>%
  subset(VWC_20cm > 0) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

P9_2021 <- P9_2021 %>%
  filter(LocalDateTime > '5/1/2021') %>%
  filter(LocalDateTime < '9/30/2021') %>%
  subset(VWC_20cm > 0) %>%
  na.omit() %>%
  select(LocalDateTime, VWC_20cm)

##join GIRF 2021 SM
GIRF2021_list <- list(P1_2021, P2_2021, P3_2021, P4_2021, P5_2021, P6_2021, P7_2021, P8_2021, P9_2021)
GIRF2021_SM_sensor_seasonal <- GIRF2021_list %>%
  reduce(full_join, by = 'LocalDateTime') %>%
  rowwise() %>%
  mutate(average = mean(c_across(starts_with("VWC")), na.rm=TRUE)) %>%
  ungroup() %>%
  get_summary_stats(average, type = 'full') 

GIRF2021_SM_sensor_oneweekbeforepulse <- GIRF2021_list %>%
  reduce(full_join, by = 'LocalDateTime') %>%
  mutate(LocalDateTime = as.Date(as.character(LocalDateTime), format = "%m/%d/%Y")) %>%
  filter(LocalDateTime > '2021/06/14' & LocalDateTime < '2021/06/20') %>%
  rowwise() %>%
  mutate(average = mean(c_across(starts_with("VWC")), na.rm=TRUE)) %>%
  ungroup() %>%
  get_summary_stats(average, type = 'full') 

GIRF2021_SM_sensor_pre_pulse <- GIRF2021_list %>%
  reduce(full_join, by = 'LocalDateTime') %>%
  mutate(LocalDateTime = as.Date(as.character(LocalDateTime), format = "%m/%d/%Y")) %>%
  filter(LocalDateTime == '2021/06/20') %>%
  rowwise() %>%
  mutate(average = mean(c_across(starts_with("VWC")), na.rm=TRUE)) %>%
  ungroup() %>%
  get_summary_stats(average, type = 'full') 

GIRF2021_SM_sensor_pulse <- GIRF2021_list %>%
  reduce(full_join, by = 'LocalDateTime') %>%
  mutate(LocalDateTime = as.Date(as.character(LocalDateTime), format = "%m/%d/%Y")) %>%
  filter(LocalDateTime > '2021/06/21' & LocalDateTime < '2021/07/08') %>%
  rowwise() %>%
  mutate(average = mean(c_across(starts_with("VWC")), na.rm=TRUE)) %>%
  ungroup() %>%
  get_summary_stats(average, type = 'full') 


##climate
TM <- read.csv("Data/RB_TM_C_SourceID_1_QC_0_Year_2021.csv")
TM_soil <- read.csv("Data/RB_TM_C_VWC_10cm_Avg_SourceID_1_QC_1.csv") 
GIRF_temp <- read.csv("Data/RB_GIRF_C_AirTemp_ST110_Avg_SourceID_1_QC_1 (1).csv") 
GIRF_precip <- read.csv("Data/RB_GIRF_C_Precip_Tot_Avg_SourceID_1_QC_1 (1).csv") 

##TM
TM$LocalDateTime <- mdy_hm(TM$LocalDateTime)
TM_soil$LocalDateTime <- mdy_hm(TM_soil$LocalDateTime)
TM_soil <- TM_soil %>%
  subset(VWC_10cm_Avg > 0) %>%
  na.omit()


TM_new <- TM[,c(1,7,40)]
TM_new <- subset(TM_new, LocalDateTime > "2021-05-01 00:00:00" & LocalDateTime < "2021-09-30 23:00:00")
TM_new$Precip_Event <- TM_new$Precip_Tot_Avg - lag(TM_new$Precip_Tot_Avg) ##find precip change
TM_new$Precip_Event <- replace(TM_new$Precip_Event, which(TM_new$Precip_Event<0), 0) ##replace negative values
TM_new <- TM_new[-1,]

TM_new1 <- merge(TM_new, TM_soil, by = "LocalDateTime")
TM_new1 <- TM_new1[,c(1,2,4,7)]


##GIRF
GIRF_precip$Precip_Event <- GIRF_precip$Precip_Tot_Avg - lag(GIRF_precip$Precip_Tot_Avg) ##find precip change
GIRF_precip$Precip_Event <- replace(GIRF_precip$Precip_Event, which(GIRF_precip$Precip_Event<0), 0) ##replace negative values
GIRF_precip <- GIRF_precip[-1, c(1,7)]
GIRF_temp <- GIRF_temp[, c(1,4)]

GIRF <-merge(GIRF_temp, GIRF_precip, by = "LocalDateTime")
GIRF$LocalDateTime <- mdy_hm(GIRF$LocalDateTime)

GIRF1 <- GIRF %>%
  subset(LocalDateTime > "2021-05-01 00:00:00" & LocalDateTime < "2021-09-30 23:00:00") %>%
  subset(AirTemp_ST110_Avg >0)%>%
  subset(Precip_Event < 100)

GIRF2 <- GIRF %>%
  subset(LocalDateTime > "2020-05-01 00:00:00" & LocalDateTime < "2020-09-30 23:00:00") %>%
  subset(AirTemp_ST110_Avg >0) %>%
  subset(Precip_Event < 100)


TM_new1 <- TM_new1 %>%
  mutate(month = month(LocalDateTime))

##seasonal
TM_new1_avg_temp <- TM_new1 %>%
  #group_by(month) %>%
  get_summary_stats(AirTemp_ST110_Avg, type = 'full') 

##pre-pulse
TM_pre_pulse_avg_temp <- TM_new1 %>%
  filter(LocalDateTime > '2021/06/17 00:00:00' & LocalDateTime < '2021/06/17 23:45:00') %>%
  get_summary_stats(AirTemp_ST110_Avg, type = 'full') 

##pulse
TM_pulse_avg_temp <- TM_new1 %>%
  filter(LocalDateTime > '2021/06/17 00:00:00' & LocalDateTime < '2021/06/27 23:45:00') %>%
  get_summary_stats(AirTemp_ST110_Avg, type = 'full') 

##average daily
TM_precip <- TM_new1 %>%
  get_summary_stats(Precip_Event, type = 'full') 
TM_new1 %>% top_n(1, Precip_Event) ##8/1/2021 19:00 had the max precip_event

##pulse total precip
TM_pulse_total_precip <- TM_new1 %>%
  filter(LocalDateTime > '2021/06/17 00:00:00' & LocalDateTime < '2021/06/27 23:45:00') %>%
  pull(Precip_Event) %>%
  sum

##total season precip
TM_prep_total <- sum(TM_new1$Precip_Event) ##total precip

##average monthly
TM_precip_monthly <- sum(TM_new1$Precip_Event)/5

##pre-pulse
TM_SM <- TM_new1 %>% 
  filter(LocalDateTime > '2021/06/23 00:00:00' & LocalDateTime < '2021/06/23 23:45:00') %>%
  get_summary_stats(VWC_10cm_Avg, type = 'full') 

##pulse
TM_SM_pre_pulse <- TM_new1 %>% 
  filter(LocalDateTime > '2021/06/17 00:00:00' & LocalDateTime < '2021/06/27 23:45:00') %>%
  get_summary_stats(VWC_10cm_Avg, type = 'full') 

##pulse
TM_SM <- TM_new1 %>%
  get_summary_stats(VWC_10cm_Avg, type = 'full') 


##GIRF 2020
GIRF2020 <- GIRF2 %>%
  mutate(month = month(LocalDateTime))

##pulse temp
GIRF2020_pulse_temp <- GIRF2020 %>%
  filter(LocalDateTime > '2020-09-03 09:00:00' & LocalDateTime < '2020-09-22 23:45:00') %>%
  get_summary_stats(AirTemp_ST110_Avg, type = 'full')

##pre-pulse temp
GIRF2020_pre_pulse_temp <- GIRF2020 %>%
  filter(LocalDateTime > '2020-09-03 00:00:00' & LocalDateTime < '2020-09-03 07:00:00') %>%
  get_summary_stats(AirTemp_ST110_Avg, type = 'full')

##seasonal temp
GIRF2020_season_temp <- GIRF2020 %>%
  get_summary_stats(AirTemp_ST110_Avg, type = 'full')

##pulse total precip
GIRF_pulse_total_precip <- GIRF2020 %>%
  filter(LocalDateTime > '2020-09-03 09:00:00' & LocalDateTime < '2020-09-22 23:45:00') %>%
  pull(Precip_Event) %>%
  sum

##prior storm events total precip
GIRF_pulse_total_precip <- GIRF2020 %>%
  filter(LocalDateTime > '2020-06-27 00:00:00' & LocalDateTime < '2020-06-30 23:45:00') %>%
  pull(Precip_Event) %>%
  sum


#daily average
GIRF2020_daily_avg <- GIRF2020 %>%
  get_summary_stats(Precip_Event, type = 'full')


##total precip
GIRF2020_total_precip <- sum(GIRF2020$Precip_Event)



#GIRF 2021
GIRF2021 <- GIRF1 %>%
  mutate(month = month(LocalDateTime))

##pulse temp
GIRF2021_pulse_temp <- GIRF2021 %>%
  filter(LocalDateTime > '2021-06-21 09:00:00' & LocalDateTime < '2021-07-07 23:45:00') %>%
  get_summary_stats(AirTemp_ST110_Avg, type = 'full')

##pre-pulse temp
GIRF2021_pre_pulse_temp <- GIRF2021 %>%
  filter(LocalDateTime > '2021-06-21 00:00:00' & LocalDateTime < '2021-06-21 07:00:00') %>%
  get_summary_stats(AirTemp_ST110_Avg, type = 'full')

##seasonal temp
GIRF2021_seasonal_temp <- GIRF2021 %>%
  get_summary_stats(AirTemp_ST110_Avg, type = 'full')

##pulse total precip
GIRF21_pulse_total_precip <- GIRF2021 %>%
  filter(LocalDateTime > '2021-06-21 09:00:00' & LocalDateTime < '2021-07-07 23:45:00') %>%
  pull(Precip_Event) %>%
  sum

##daily average
GIRF2021_1 <- GIRF2021 %>%
  get_summary_stats(Precip_Event, type = 'full')

##average monthly
GIRF2021_monthly_precip <- sum(GIRF2021$Precip_Event)/5

##total precip
GIRF2021_total_precip <- sum(GIRF2021$Precip_Event)


