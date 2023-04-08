##GIRF 2020-2021 and TM 2021 PSEM
##Authored by: Yvette  D. Hastings
##Written: May 25, 2022
##Updated: July 14, 2022

library(readxl)
library(piecewiseSEM)
library(dplyr)
library(tidyverse)
library(lavaan)
library(lavaanPlot)
library(semEff)
library(corrplot)


all_pulse <- read_excel("2020-2021 Soils Master File.xlsx", sheet = 'Pulse Data Masterfile (R)')

##June 2021
June2021 <- filter(all_pulse, sampling.campaign == "June 2021 GIRF Pulse")  %>%
  drop_na(mb.doc)

June2021$date <- ifelse(June2021$sampling.date == as.Date('2021-06-17'), 1,
                            ifelse(June2021$sampling.date == as.Date('2021-06-21'), 2,
                                   ifelse(June2021$sampling.date == as.Date('2021-06-22'), 3,
                                          ifelse(June2021$sampling.date == as.Date('2021-06-23'), 4,
                                                 ifelse(June2021$sampling.date == as.Date('2021-06-24'),5,
                                                        ifelse(June2021$sampling.date == as.Date('2021-06-28'),6,
                                                               ifelse(June2021$sampling.date == as.Date('2021-06-30'),7,
                                                                      ifelse(June2021$sampling.date == as.Date('2021-06-25'),8,
                                                                             ifelse(June2021$sampling.date == as.Date('2021-06-27'),9,10)))))))))
June2021$plot.type <- ifelse(June2021$treatment == 'Mulch', 1,
                                 ifelse(June2021$treatment =='Grass',2,
                                        ifelse(June2021$treatment =='Diverse',3,4)))

June2021_1 <- subset(June2021, select = -c(1:3, 5:7, 13:14, 23:24, 27))
colnames(June2021_1) <- c("Plot Number", "Bulk Density", "pH", "Total N", "Inorganic N", "Organic N", "Gravimetric Soil Moisture", "Organic Matter Content", 
                        "Microbial Biomass C", "Microbial Biomass N", "LAP", "AP", "BG", "POX", "Proteolytic Rate", "Protein Turnover Rate", 
                        "Native Protein Concentration" ,"Date", "Plot Type")

##Sept 2020
Sept2020 <- filter(all_pulse, sampling.campaign == "GIRF September 2020 Pulse")

Sept2020$date <- ifelse(Sept2020$sampling.date == as.Date('2020-09-03'), 1,
                            ifelse(Sept2020$sampling.date == as.Date('2020-09-09'), 2,
                                   ifelse(Sept2020$sampling.date == as.Date('2020-09-11'), 3,
                                          ifelse(Sept2020$sampling.date == as.Date('2020-09-15'), 4,5))))
Sept2020$plot.type <- ifelse(Sept2020$treatment == 'Mulch', 1,
                                 ifelse(Sept2020$treatment =='Grass',2,3))

Sept2020_1 <- subset(Sept2020, select = -c(1:3, 5:7, 13:14, 23:24, 27))

colnames(Sept2020_1) <- c("Plot Number", "Bulk Density", "pH", "Total N", "Inorganic N", "Organic N", "Gravimetric Soil Moisture", "Organic Matter Content", 
                        "Microbial Biomass C", "Microbial Biomass N", "LAP", "AP", "BG", "POX", "Proteolytic Rate", "Protein Turnover Rate", 
                        "Native Protein Concentration" ,"Date", "Plot Type")

##TM2021
TM2021 <- filter(all_pulse, sampling.campaign == "June 2021 TM Natural Pulse")

TM2021$date <- ifelse(TM2021$sampling.date == as.Date('2021-06-17'), 1,
                          ifelse(TM2021$sampling.date == as.Date('2021-06-25'), 2, 3))
TM2021$treatment[TM2021$treatment == 'Reference'] <- 1

TM2021_1 <- subset(TM2021, select = -c(1:3, 5:8, 13:14, 23:24, 27))

colnames(TM2021_1) <- c("Plot Number", "pH", "Total N", "Inorganic N", "Organic N", "Gravimetric Soil Moisture", "Organic Matter Content", 
                      "Microbial Biomass C", "Microbial Biomass N", "LAP", "AP", "BG", "POX", "Proteolytic Rate", "Protein Turnover Rate", 
                      "Native Protein Concentration" ,"Date")


##Sept 2020 Correlation matrix
# Sept2020_COR[,chars] <- as.data.frame(apply(Sept2020_COR[,chars], 2, as.numeric))

##correlation matrix
S <- cor(Sept2020_1)

##assign p-values to correlation
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
S_p.mat <- cor.mtest(Sept2020_1, method = 'pearson')

##correlation plot with significant only
sept <- corrplot.mixed(S, upper = "ellipse", lower = "number", tl.pos = "lt",
         p.mat = S_p.mat, sig.level = 0.05, insig = "blank", tl.col = 'black',
         order = 'alphabet', addgrid.col = 'grey', 
         mar=c(0,0,3,0))

write.csv(S, "SeptCorr.csv")
##PSEM
##fit models - Sept 2020
sept.om <- lm(percent.organic.matter ~ gravimetric.sm + enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = Sept2020)
sept.doc <- lm(mb.doc ~ gravimetric.sm, data = Sept2020)
sept.tdn <- lm(mb.tdn ~ gravimetric.sm, data = Sept2020)

sept.lap <- lm(enzyme.lap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = Sept2020)
sept.ap <- lm(enzyme.ap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = Sept2020)
sept.bg <- lm(enzyme.bg ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = Sept2020)
sept.pox <- lm(enzyme.pox ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = Sept2020)

sept.pr <- lm(potential.net.proteolytic.rate ~ percent.organic.matter + enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = Sept2020)

sept.on <- lm(organic.nitrogen ~ potential.net.proteolytic.rate + percent.organic.matter + 
              enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = Sept2020)
sept.in <- lm(inorganic.nitrogen ~ organic.nitrogen + enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = Sept2020)

sept.ph <- lm(ph ~ inorganic.nitrogen, data = Sept2020)

Sept_2020.psem <- psem(## Regressions
  sept.om,
  sept.doc, 
  sept.tdn, 
  sept.lap, 
  sept.ap, 
  sept.bg, 
  sept.pox, 
  sept.pr, 
  sept.on, 
  sept.in, 
  sept.ph,
  ## Covariances
  ph %~~% gravimetric.sm,
  percent.organic.matter %~~% mb.doc,
  percent.organic.matter %~~% mb.tdn,
  mb.doc %~~% percent.organic.matter,
  mb.tdn %~~% percent.organic.matter,
  mb.doc %~~% mb.tdn,
  mb.tdn %~~% mb.doc,
  mb.doc %~~% inorganic.nitrogen,
  mb.tdn %~~% inorganic.nitrogen,
  enzyme.lap %~~% enzyme.ap,
  enzyme.lap %~~% enzyme.bg,
  enzyme.lap %~~% enzyme.pox,
  enzyme.ap %~~% enzyme.bg,
  enzyme.ap %~~% enzyme.pox,
  enzyme.ap %~~% enzyme.lap,
  enzyme.bg %~~% enzyme.ap,
  enzyme.bg %~~% enzyme.lap,
  enzyme.bg %~~% enzyme.pox,
  enzyme.pox %~~% enzyme.ap,
  enzyme.pox %~~% enzyme.lap,
  enzyme.pox %~~% enzyme.bg,
  ## Data
  data = as.data.frame(Sept2020))

summary(Sept_2020.psem)

##test to see which inaction of BG and AP influence one of the linking variables
summary(psem(lm(percent.organic.matter ~ enzyme.ap*enzyme.bg, data = Sept2020)))

# sem_graph <- plot(Sept_2020.psem, return = TRUE)
# DiagrammeR::render_graph(sem_graph)

## To change the layout, print out the node table
sem_graph$nodes_df

## Updated node positions:
sem_graph <- plot(Sept_2020.psem, digits = 3, node_attrs = list(
  shape = "ellipse", color = "black", 
  height = 0.5, width = 1, cex = 0.5,
  fillcolor = "grey",
  ##  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)
  x = c( 2, 2, 2, 3, 3, 3, 3, 5, 6, 8.5, 0, 11, 11, 10, 12, 10, 9, 10, 11, 0, 9, 12),
  y = c( 3, 5, 7, 2, 4, 6, 8, 3, 5, 5, 4, 7, 6, 5, 5, 4, 3, 2, 3, 6, 7, 4)
), 
return = TRUE)
DiagrammeR::render_graph(sem_graph)

##calculate effects
##bootstrap and save standarized direct effects
##https://murphymv.github.io/semEff/articles/semEff.html

sept.sem <- list(
  lm(percent.organic.matter ~ gravimetric.sm + enzyme.lap + enzyme.ap * enzyme.bg + enzyme.pox, data = Sept2020),
  lm(mb.doc ~ gravimetric.sm , data = Sept2020),
  lm(mb.tdn ~ gravimetric.sm , data = Sept2020),
  lm(mb.doc ~ inorganic.nitrogen, data = Sept2020),
  lm(mb.tdn ~ inorganic.nitrogen, data = Sept2020),
  lm(enzyme.lap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = Sept2020),
  lm(enzyme.ap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = Sept2020),
  lm(enzyme.bg ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = Sept2020),
  lm(enzyme.pox ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = Sept2020),
  lm(potential.net.proteolytic.rate ~ percent.organic.matter + enzyme.lap + enzyme.ap * enzyme.bg + enzyme.pox, data = Sept2020),
  lm(organic.nitrogen ~ potential.net.proteolytic.rate + percent.organic.matter +
       enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = Sept2020),
  lm(inorganic.nitrogen ~ organic.nitrogen + enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = Sept2020)
)

system.time(
  sept.sem.boot <- bootEff(sept.sem, R = 100, seed = 13, parallel = "no")
  
)

##caluculate effects
(sept.sem.eff <- semEff(sept.sem.boot))
summary(sept.sem.eff)
summary(sept.sem.eff, response = "mb.doc")
summary(sept.sem.eff, response = "mb.tdn")
summary(sept.sem.eff, response = "enzyme.lap")
summary(sept.sem.eff, response = "enzyme.ap")
summary(sept.sem.eff, response = "enzyme.bg")
summary(sept.sem.eff, response = "enzyme.pox")
summary(sept.sem.eff, response = "percent.organic.matter")

##find path mediators for indirect path predictor, update predictor and response variables to find mediator effects
summary(
  semEff(sept.sem.boot, predictor = "gravimetric.sm"),
  response = "enzyme.lap"
)


##June 2021 Correlation matrix
##correlation matrix
J <- cor(June2021_1)

##assign p-values to correlation
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
J_p.mat <- cor.mtest(June2021_1, method = 'pearson')

##correlation plot with significant only
june <- corrplot.mixed(J, upper = "ellipse", lower = "number", tl.pos = "lt",
                       p.mat = J_p.mat, sig.level = 0.05, insig = "blank", tl.col = 'black',
                       order = 'alphabet', addgrid.col = 'grey', 
                       mar=c(0,0,3,0))

write.csv(J, "JuneCorr.csv")

##PSEM
##fit models - June 2021
june.om <- lm(percent.organic.matter ~ gravimetric.sm + enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = June2021)
june.doc <- lm(mb.doc ~ gravimetric.sm, data = June2021)
june.tdn <- lm(mb.tdn ~ gravimetric.sm, data = June2021)

june.lap <- lm(enzyme.lap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = June2021)
june.ap <- lm(enzyme.ap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = June2021)
june.bg <- lm(enzyme.bg ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = June2021)
june.pox <- lm(enzyme.pox ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = June2021)

june.pr <- lm(potential.net.proteolytic.rate ~ percent.organic.matter + 
                enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = June2021)

june.on <- lm(organic.nitrogen ~ potential.net.proteolytic.rate + percent.organic.matter + 
                enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = June2021)
june.in <- lm(inorganic.nitrogen ~ organic.nitrogen + enzyme.lap + enzyme.ap + 
                enzyme.bg + enzyme.pox, data = June2021)

june.ph <- lm(ph ~ inorganic.nitrogen, data = June2021)


June_2021.psem <- psem(## Regressions
  june.om,
  june.doc, 
  june.tdn, 
  june.lap, 
  june.ap, 
  june.bg, 
  june.pox, 
  june.pr, 
  june.on, 
  june.in, 
  june.ph,
  ## Covariances
  ph %~~% gravimetric.sm,
  percent.organic.matter %~~% mb.doc,
  percent.organic.matter %~~% mb.tdn,
  mb.doc %~~% percent.organic.matter,
  mb.tdn %~~% percent.organic.matter,
  mb.doc %~~% mb.tdn,
  mb.tdn %~~% mb.doc,
  mb.doc %~~% inorganic.nitrogen,
  mb.tdn %~~% inorganic.nitrogen,
  enzyme.lap %~~% enzyme.ap,
  enzyme.lap %~~% enzyme.bg,
  enzyme.lap %~~% enzyme.pox,
  enzyme.ap %~~% enzyme.bg,
  enzyme.ap %~~% enzyme.pox,
  enzyme.ap %~~% enzyme.lap,
  enzyme.bg %~~% enzyme.ap,
  enzyme.bg %~~% enzyme.lap,
  enzyme.bg %~~% enzyme.pox,
  enzyme.pox %~~% enzyme.ap,
  enzyme.pox %~~% enzyme.lap,
  enzyme.pox %~~% enzyme.bg,
  ## Data
  data = as.data.frame(June2021))

summary(June_2021.psem)


## To change the layout, print out the node table
sem_graph$nodes_df

## Updated node positions:
sem_graph <- plot(June_2021.psem, digits = 3, node_attrs = list(
  shape = "ellipse", color = "black", 
  height = 0.5, width = 1, cex = 0.5,
  fillcolor = "orange",
  ##  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22
  x = c( 2, 2, 2, 3, 3, 3, 3, 5, 6, 8.5, 0, 11, 11, 10, 12, 10, 9, 10, 11, 0, 9, 12),
  y = c( 3, 5, 7, 2, 4, 6, 8, 3, 5, 5, 4, 7, 6, 5, 5, 4, 3, 2, 3, 6, 7, 4)
), 
return = TRUE)
DiagrammeR::render_graph(sem_graph)


##calculate effects
##bootstrap and save standarized direct effects
##https://murphymv.github.io/semEff/articles/semEff.html

june.sem <- list(
  lm(percent.organic.matter ~ gravimetric.sm, data = June2021),
  lm(mb.doc ~ gravimetric.sm, data = June2021),
  lm(mb.tdn ~ gravimetric.sm, data = June2021),
  lm(enzyme.lap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = June2021),
  lm(enzyme.ap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = June2021),
  lm(enzyme.bg ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = June2021),
  lm(enzyme.pox ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = June2021),
  lm(potential.net.proteolytic.rate ~ percent.organic.matter + 
       enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = June2021),
  lm(organic.nitrogen ~ potential.net.proteolytic.rate + percent.organic.matter +
       enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = June2021),
  lm(inorganic.nitrogen ~ organic.nitrogen + enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox +ph, data = June2021)
)

system.time(
  june.sem.boot <- bootEff(june.sem, R = 100, seed = 13, parallel = "no")
  
)

##caluculate effects
(june.sem.eff <- semEff(june.sem.boot))
summary(june.sem.eff)
summary(june.sem.eff, response = "inorganic.nitrogen")
summary(june.sem.eff, response = "mb.doc")
summary(june.sem.eff, response = "mb.tdn")
summary(june.sem.eff, response = "enzyme.lap")
summary(june.sem.eff, response = "enzyme.ap")
summary(june.sem.eff, response = "enzyme.bg")
summary(june.sem.eff, response = "enzyme.pox")

##find path mediators for indirect path predictor
summary(
  semEff(june.sem.boot, predictor = "ph"),
  response = "inorganic.nitrogen"
)

summary(
  semEff(june.sem.boot, predictor = "enzyme.lap"),
  response = "inorganic.nitrogen"
)

summary(
  semEff(june.sem.boot, predictor = "enzyme.ap"),
  response = "inorganic.nitrogen"
)

summary(
  semEff(june.sem.boot, predictor = "mb.tdn"),
  response = "inorganic.nitrogen"
)


##TM 2021 Correlation matrix
# Sept2020_COR[,chars] <- as.data.frame(apply(Sept2020_COR[,chars], 2, as.numeric))

##correlation matrix
TM <- cor(TM2021_1)

##assign p-values to correlation
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
TM_p.mat <- cor.mtest(TM2021_1, method = 'pearson') ##remove bulk density and volumetric.sm since observations are missing

##correlation plot with significant only
tm <- corrplot.mixed(TM, upper = "ellipse", lower = "number", tl.pos = "lt",
                     p.mat = TM_p.mat, sig.level = 0.05, insig = "blank", tl.col = 'black',
                     order = 'alphabet', addgrid.col = 'grey', 
                     mar=c(0,0,3,0))

write.csv(TM, "TMCorr.csv")

##PSEM

##fit models - TM 2021
tm.om <- lm(percent.organic.matter ~ gravimetric.sm + enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = TM2021)
tm.doc <- lm(mb.doc ~ gravimetric.sm, data = TM2021)
tm.tdn <- lm(mb.tdn ~ gravimetric.sm, data = TM2021)

tm.lap <- lm(enzyme.lap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = TM2021)
tm.ap <- lm(enzyme.ap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = TM2021)
tm.bg <- lm(enzyme.bg ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = TM2021)
tm.pox <- lm(enzyme.pox ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = TM2021)

tm.pr <- lm(potential.net.proteolytic.rate ~ percent.organic.matter + 
                enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = TM2021)

tm.on <- lm(organic.nitrogen ~ potential.net.proteolytic.rate + percent.organic.matter + 
                enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = TM2021)
tm.in <- lm(inorganic.nitrogen ~ organic.nitrogen + enzyme.lap + enzyme.ap + 
                enzyme.bg + enzyme.pox, data = TM2021)

tm.ph <- lm(ph ~ inorganic.nitrogen, data = June2021)


TM_2021.psem <- psem(## Regressions
  tm.om,
  tm.doc, 
  tm.tdn, 
  tm.lap, 
  tm.ap, 
  tm.bg, 
  tm.pox, 
  tm.pr, 
  tm.on, 
  tm.in, 
  tm.ph,
  ## Covariances
  ph %~~% gravimetric.sm,
  percent.organic.matter %~~% mb.doc,
  percent.organic.matter %~~% mb.tdn,
  mb.doc %~~% percent.organic.matter,
  mb.tdn %~~% percent.organic.matter,
  mb.doc %~~% mb.tdn,
  mb.tdn %~~% mb.doc,
  mb.doc %~~% inorganic.nitrogen,
  mb.tdn %~~% inorganic.nitrogen,
  enzyme.lap %~~% enzyme.ap,
  enzyme.lap %~~% enzyme.bg,
  enzyme.lap %~~% enzyme.pox,
  enzyme.ap %~~% enzyme.bg,
  enzyme.ap %~~% enzyme.pox,
  enzyme.ap %~~% enzyme.lap,
  enzyme.bg %~~% enzyme.ap,
  enzyme.bg %~~% enzyme.lap,
  enzyme.bg %~~% enzyme.pox,
  enzyme.pox %~~% enzyme.ap,
  enzyme.pox %~~% enzyme.lap,
  enzyme.pox %~~% enzyme.bg,
  ## Data
  data = as.data.frame(TM2021))

summary(TM_2021.psem) ##not working, run summary on individual models

summary(psem(lm(enzyme.ap ~ enzyme.bg, data = TM2021))) ##lm models on those as as covariances in psem model give a slightly differet result. Use value from actual model but use this to get significance.

## To change the layout, print out the node table
sem_graph$nodes_df

## Updated node positions:
sem_graph <- plot(TM_2021.psem, digits = 3, node_attrs = list(
  shape = "ellipse", color = "black", 
  height = 0.5, width = 1, cex = 0.5,
  fillcolor = "blue", 
  ##  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22
  x = c( 2, 2, 2, 3, 3, 3, 3, 5, 6, 8.5, 0, 11, 11, 10, 12, 10, 9, 10, 11, 0, 9, 12),
  y = c( 3, 5, 7, 2, 4, 6, 8, 3, 5, 5, 4, 7, 6, 5, 5, 4, 3, 2, 3, 6, 7, 4)
), 
return = TRUE)
DiagrammeR::render_graph(sem_graph)


##calculate effects
##bootstrap and save standarized direct effects
##https://murphymv.github.io/semEff/articles/semEff.html

tm.sem <- list(
  lm(percent.organic.matter ~ gravimetric.sm, data = TM2021),
  lm(mb.doc ~ gravimetric.sm, data = TM2021),
  lm(mb.tdn ~ gravimetric.sm, data = TM2021),
  lm(enzyme.lap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = TM2021),
  lm(enzyme.ap ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = TM2021),
  lm(enzyme.bg ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = TM2021),
  lm(enzyme.pox ~ gravimetric.sm + ph + mb.doc + mb.tdn, data = TM2021),
  lm(potential.net.proteolytic.rate ~ percent.organic.matter + 
       enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = TM2021),
  lm(organic.nitrogen ~ potential.net.proteolytic.rate + percent.organic.matter +
       enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = TM2021),
  lm(inorganic.nitrogen ~ organic.nitrogen + enzyme.lap + enzyme.ap + enzyme.bg + enzyme.pox, data = TM2021)
)

system.time(
  tm.sem.boot <- bootEff(tm.sem, R = 100, seed = 13, parallel = "no")
  
)

##caluculate effects
(tm.sem.eff <- semEff(tm.sem.boot))
summary(tm.sem.eff)
summary(tm.sem.eff, response = "inorganic.nitrogen")
summary(tm.sem.eff, response = "mb.doc")
summary(tm.sem.eff, response = "mb.tdn")
summary(tm.sem.eff, response = "enzyme.lap")
summary(tm.sem.eff, response = "enzyme.ap")
summary(tm.sem.eff, response = "enzyme.bg")
summary(tm.sem.eff, response = "enzyme.pox")

##find path mediators for indirect path predictor
summary(
  semEff(tm.sem.boot, predictor = "gravimetric.sm"),
  response = "inorganic.nitrogen"
)




