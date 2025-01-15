### load the libraries
library(lme4)
library(emmeans)
library(sjPlot)
library(report)
library(ggplot2)
foldN = "."
wtMice <- c(11,18,21,24,45,46,47)
pvGSMice <- c(57,61,62,63,65,4)
pvFamStimMice <- c(48,50,52,53,54,57,61,62,63,65,4)
######################## ripple properties #####################################
rm()
dataStrings <- c("SWR_rate","SWR_duration","SWR_power","SWR_coactivProb")

for (dt in dataStrings) {
  
  tblFilename = paste0(foldN, "/TableData_", dt, ".txt")
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create .txt file to save LMM outputs 
  
  mydata <- read.csv(tblFilename, head=T)
  mydata <- mydata[!is.na(mydata$DataY), ]
  
  #factorize variables to be treated as categories/levels instead of numerical value
  mydata$AnimalID = as.factor(mydata$AnimalID)
  mydata$NovelDay = as.factor(mydata$NovelDay)
  mydata$Environ = as.factor(mydata$Environ)
  mydata$StimCond = as.factor(mydata$StimCond)
  
  sink(txtFilename)
  wt <- mydata[is.element(mydata$AnimalID, wtMice), ] 
  obj.lmer = lmerTest::lmer(DataY ~ NovelDay*Environ + (1|AnimalID), data = wt)
  print(summary(obj.lmer))
  print(anova(obj.lmer))
  print(contrast(emmeans(obj.lmer, specs=c("NovelDay","Environ")), "pairwise", adjust = "tukey"))
  
  wt <- mydata[is.element(mydata$AnimalID, wtMice), ] 
  obj.lmer = lmerTest::lmer(DataY ~ Environ + (1|AnimalID), data = wt)
  print(summary(obj.lmer))
  print(anova(obj.lmer))
  print(contrast(emmeans(obj.lmer, specs=c("Environ")), "pairwise", adjust = "tukey"))
  
  pv_nov <- mydata[is.element(mydata$AnimalID, pvGSMice) & mydata$Environ == 2, ] 
  obj.lmer = lmerTest::lmer(DataY ~ NovelDay*StimCond + (1|AnimalID), data = pv_nov)
  print(summary(obj.lmer))
  print(anova(obj.lmer))
  print(contrast(emmeans(obj.lmer, specs=c("NovelDay","StimCond")), "pairwise", adjust = "tukey"))
  
  pv_nov <- mydata[is.element(mydata$AnimalID, pvGSMice) & mydata$Environ == 2, ] 
  obj.lmer = lmerTest::lmer(DataY ~ StimCond + (1|AnimalID), data = pv_nov)
  print(summary(obj.lmer))
  print(anova(obj.lmer))
  print(contrast(emmeans(obj.lmer, specs=c("StimCond")), "pairwise", adjust = "tukey"))
  
  pv <- mydata[is.element(mydata$AnimalID, pvGSMice) & mydata$Environ != 2, ] 
  obj.lmer = lmerTest::lmer(DataY ~ NovelDay*StimCond + (1|AnimalID), data = pv)
  print(summary(obj.lmer))
  print(anova(obj.lmer))
  print(contrast(emmeans(obj.lmer, specs=c("NovelDay","StimCond")), "pairwise", adjust = "tukey"))
  
  pv <- mydata[is.element(mydata$AnimalID, pvGSMice) & mydata$Environ != 2, ] 
  obj.lmer = lmerTest::lmer(DataY ~ StimCond + (1|AnimalID), data = pv)
  print(summary(obj.lmer))
  print(anova(obj.lmer))
  print(contrast(emmeans(obj.lmer, specs=c("StimCond")), "pairwise", adjust = "tukey"))
  
  sink()
}


