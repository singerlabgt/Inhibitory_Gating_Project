library(lme4)
library(emmeans)
library(sjPlot)
library(report)
library(ggplot2)
library(pbkrtest)
setwd("ExtendedDataFigure03")

wtMice <- c(11,18,21,24,45,46,47)
pvGSMice <- c(57,61,62,63,65)
pvFamStimMice <- c(48,50,52,53,54,57,61,62,63,65)
dataStrings <- c("FRxSpeedQuartile")
cellTypes <-c("NarrowInt","WideInt","Pyr") #corresponds to 1-3

for (dt in dataStrings)  {
  
  
  #load MATLAB output (in table format in .txt) file in R
  datadir <- "."
  foldN = "."
  
  tblFilename = paste0(datadir, "/TableData_", dt, ".txt")
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create .txt file to save LMM outputs 
  
  mydata <- read.csv(tblFilename, head=T)
  mydata <- mydata[!is.na(mydata$CellType) & !is.na(mydata$DataY),]
  
  #factorize variables to be treated as categories/levels instead of numerical value
  mydata$AnimalID = as.factor(mydata$AnimalID)
  mydata$SpeedQuartile = as.factor(mydata$SpeedQuartile)
  mydata$RecDate = as.factor(mydata$RecDate)
  mydata$NovelDay = as.factor(mydata$NovelDay)
  mydata$CellType = as.factor(mydata$CellType)
  mydata$CellID = as.factor(mydata$CellID)
  
  #create model
  wtFRxSpeedQ_NSint<- mydata[is.element(mydata$AnimalID, wtMice) & mydata$CellType == 1, ] #NS int from WT mice only
  obj.lmer1 = lmerTest::lmer(DataY ~ SpeedQuartile + (1|AnimalID/CellID), data = wtFRxSpeedQ_NSint)
  
  wtFRxSpeedQ_WSint<- mydata[is.element(mydata$AnimalID, wtMice) & mydata$CellType == 2, ] #WS int from WT mice only
  obj.lmer2 = lmerTest::lmer(DataY ~ SpeedQuartile + (1|AnimalID/CellID), data = wtFRxSpeedQ_WSint)
  
  wtFRxSpeedQ_pyr<- mydata[is.element(mydata$AnimalID, wtMice) & mydata$CellType == 3, ] #WS int from WT mice only
  obj.lmer3 = lmerTest::lmer(DataY ~ SpeedQuartile + (1|AnimalID/CellID), data = wtFRxSpeedQ_pyr)
  
  wtFRxSpeedQ_all<- mydata[is.element(mydata$AnimalID, wtMice) , ] #WS int from WT mice only
  obj.lmer4 = lmerTest::lmer(DataY ~ SpeedQuartile*CellType + (1|AnimalID/CellID), data = wtFRxSpeedQ_all)
  
  ## output all stats details
  sink(txtFilename)
  
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer1, specs=c("SpeedQuartile")), "pairwise", adjust = "tukey"))
  
  print(summary(obj.lmer2))
  print(anova(obj.lmer2, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer2, specs=c("SpeedQuartile")), "pairwise", adjust = "tukey"))
  
  print(summary(obj.lmer3))
  print(anova(obj.lmer3, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer3, specs=c("SpeedQuartile")), "pairwise", adjust = "tukey"))
  
  print(summary(obj.lmer4))
  print(anova(obj.lmer4, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer4, specs=c("SpeedQuartile","CellType")), "pairwise", adjust = "tukey"))
  
  sink()
}
