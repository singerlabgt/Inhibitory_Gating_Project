library(lme4)
library(emmeans)
library(sjPlot)
library(report)
library(ggplot2)

wtMice <- c(11,18,21,24,45,46,47)
pvGSMice <- c(57,61,62,63,65,4)
pvFamStimMice <- c(48,50,52,53,54,57,61,62,63,65,4)
dataStrings <- c("spatialInfo","ratemapCorr")
cellTypes <-c("NarrowInt","WideInt","Pyr") #corresponds to 1-3

foldN = "." 
behaviorTypes = c("stimSpeed","stimLickR")
for (dt in c(1,2)) {
  
  tblFilename = paste0(foldN, "/TableData_", "normBehxIntensity", ".txt")
  txtFilename = paste0(foldN, "/stats_lme4_", behaviorTypes[dt], ".txt") #create .txt file to save LMM outputs 
  
  mydata <- read.csv(tblFilename, head=T)
  mydata <- mydata[!is.na(mydata$DataY) & mydata$BehType == dt, ]
  
  #factorize variables to be treated as categories/levels instead of numerical value
  mydata$AnimalID = as.factor(mydata$AnimalID)
  mydata$NovelDay = as.factor(mydata$NovelDay)
  mydata$Environ = as.factor(mydata$Environ)
  mydata$StimCond = as.factor(mydata$StimCond)
  mydata$StimIntensity = as.factor(mydata$StimIntensity)
  
  sink(txtFilename)
  
  pv_famGS <- mydata[is.element(mydata$AnimalID, pvGSMice) & mydata$Environ == 1 & is.element(mydata$StimCond, c(1)), ] 
  obj.lmer = lmerTest::lmer(DataY ~ StimIntensity + (1|AnimalID), data = pv_famGS)
  print(paste0(behaviorTypes[dt], " - Fam"))
  print(summary(obj.lmer))
  print(anova(obj.lmer))
  print(contrast(emmeans(obj.lmer, specs=c("StimIntensity")), "pairwise", adjust = "tukey"))
  
  for (nD in c(1,2,3)) {
    pv_nov <- mydata[is.element(mydata$AnimalID, pvGSMice) & mydata$Environ == 2 & is.element(mydata$StimCond, c(2,3)) & mydata$NovelDay == nD, ] 
    obj.lmer = lmerTest::lmer(DataY ~ StimIntensity + (1|AnimalID), data = pv_nov)
    print(paste0(behaviorTypes[dt], " - Nov", nD))
    print(summary(obj.lmer))
    print(anova(obj.lmer))
    print(contrast(emmeans(obj.lmer, specs=c("StimIntensity")), "pairwise", adjust = "tukey"))
    
  }
  
  sink()
}