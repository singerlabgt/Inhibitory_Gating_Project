### load the libraries
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

for (dt in dataStrings) {
  
  
  #load MATLAB output (in table format in .txt) file in R

  foldN = "."
  
  tblFilename = paste0(foldN, "/TableData_", dt, ".txt")
  txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create .txt file to save LMM outputs 
  
  mydata <- read.csv(tblFilename, head=T)
  mydata <- mydata[!is.na(mydata$DataY),]
  
  #factorize variables to be treated as categories/levels instead of numerical value
  mydata$AnimalID = as.factor(mydata$AnimalID)
  mydata$RecDate = as.factor(mydata$RecDate)
  mydata$NovelDay = as.factor(mydata$NovelDay)
  mydata$Environ = as.factor(mydata$Environ)
  mydata$CellID = as.factor(mydata$CellID)
  mydata$CellType = as.factor(mydata$CellType)
  mydata$StimCond = as.factor(mydata$StimCond)
  
  ## PVxAi32: compare goal stim and sham stim in novel environment only 
  #compare goal-modulated pyramidal cells only in the novel environment  
  pvGoalModPyr <- mydata[is.element(mydata$AnimalID, pvGSMice) &
                           mydata$GoalMod == 1 & #goal modulated 
                           mydata$CellType == 3 & #pyramidal cells only 
                           mydata$Environ == 2 & #novel only
                           is.element(mydata$StimCond, c(2,3)), ] #goal or sham conditions only
  obj.lmer1 = lmerTest::lmer(DataY ~ NovelDay*StimCond + (1|AnimalID/CellID), data = pvGoalModPyr)
  
  ## PVxAi32: compare non-goal modulated pyramidal cells only 
  pvNonGoalModPyr <- mydata[is.element(mydata$AnimalID, pvGSMice) &
                              mydata$GoalMod == 0 & #non-goal modulated 
                              mydata$CellType == 3 & #pyramidal cells only 
                              mydata$Environ == 2 & #novel only
                              is.element(mydata$StimCond, c(2,3)), ] 
  obj.lmer2 = lmerTest::lmer(DataY ~ NovelDay*StimCond + (1|AnimalID/CellID), data = pvNonGoalModPyr)
  
  
  ## PVxAi32: compare goal modulated narrow interneurons only 
  pvGoalModInt <- mydata[is.element(mydata$AnimalID, pvGSMice) &
                           mydata$GoalMod == 1 & #goal modulated 
                           mydata$CellType == 1 & #narrow interneurons only
                           mydata$Environ == 2 & #novel only
                           is.element(mydata$StimCond, c(2,3)), ] 
  obj.lmer3 = lmerTest::lmer(DataY ~ NovelDay*StimCond + (1|AnimalID/CellID), data = pvGoalModInt)
  
  ## PVxAi32: compare non-goal modulated pramidal cells only 
  pvNonGoalModInt <- mydata[is.element(mydata$AnimalID, pvGSMice) &
                              mydata$GoalMod == 0 & #non-goal modulated 
                              mydata$CellType == 1 & #narrow interneurons only
                              mydata$Environ == 2 & #novel only
                              is.element(mydata$StimCond, c(2,3)), ] 
  obj.lmer4 = lmerTest::lmer(DataY ~ NovelDay*StimCond + (1|AnimalID/CellID), data = pvNonGoalModInt)
  
  
  ## helpful plotting lines for future reference
  # ggplot(dat[dat$StimCond==3,],aes(x=NovelDay,y=DataY)) + geom_point() + 
  #   geom_smooth(method = "lm", level = 0.95) + geom_point() + facet_wrap(~AnimalID, nrow = 3, ncol = 2)
  
  #obj.lmer.ml=lme4::lmer(DataY ~ NovelDay*StimCond +(1|AnimalID/CellID), data=dat, REML=F) 
  #obj.lmer0.ml=lme4::lmer(DataY ~ 1+(1|AnimalID/CellID), data=dat, REML=F) 
  
  
  ################################################################################
  ################################################################################
  ## WT mice only 
  
  ## WT: compare Fam vs Nov 
  #compare goal-modulated pyramidal cells only 
  wtGoalModPyr <- mydata[is.element(mydata$AnimalID, wtMice) &
                           mydata$GoalMod == 1 & #goal modulated 
                           mydata$CellType == 3 & #pyramidal cells only
                           is.element(mydata$Environ, c(1,2)), ] 
  obj.lmer5 = lmerTest::lmer(DataY ~ NovelDay*Environ + (1|AnimalID/CellID), data = wtGoalModPyr)
  
  ## WT: compare non-goal-modulated pyramidal cells only 
  wtNonGoalModPyr <- mydata[is.element(mydata$AnimalID, wtMice) &
                              mydata$GoalMod == 0 &
                              mydata$CellType == 3 & #pyramidal cells only
                              is.element(mydata$Environ, c(1,2)), ] 
  obj.lmer6 = lmerTest::lmer(DataY ~ NovelDay*Environ + (1|AnimalID/CellID), data = wtNonGoalModPyr)
  
  ## WT: compare goal-modulated narrow interneurons only 
  wtGoalModInt <- mydata[is.element(mydata$AnimalID, wtMice) &
                           mydata$GoalMod == 1 &
                           mydata$CellType == 1 & #narrow interneurons
                           is.element(mydata$Environ, c(1,2)), ] 
  obj.lmer7 = lmerTest::lmer(DataY ~ NovelDay*Environ + (1|AnimalID/CellID), data = wtGoalModInt)
  
  ## WT: compare non-goal-modulated narrow interneurons only 
  wtNonGoalModInt <- mydata[is.element(mydata$AnimalID, wtMice) &
                              mydata$GoalMod == 0 &
                              mydata$CellType == 1 & #narrow interneurons
                              is.element(mydata$Environ, c(1,2)), ] 
  obj.lmer8 = lmerTest::lmer(DataY ~ NovelDay*Environ + (1|AnimalID/CellID), data = wtNonGoalModInt)
  
  
  
  ## output all stats details
  sink(txtFilename)
  
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer1, specs=c("NovelDay","StimCond")), "pairwise", adjust = "tukey"))
  
  print(summary(obj.lmer2))
  print(anova(obj.lmer2, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer2, specs=c("NovelDay","StimCond")), "pairwise", adjust = "tukey"))
  
  print(summary(obj.lmer3))
  print(anova(obj.lmer3, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer3, specs=c("NovelDay","StimCond")), "pairwise", adjust = "tukey"))
  
  print(summary(obj.lmer4))
  print(anova(obj.lmer4, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer4, specs=c("NovelDay","StimCond")), "pairwise", adjust = "tukey"))
  
  print(summary(obj.lmer5))
  print(anova(obj.lmer5, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer5, specs=c("NovelDay","Environ")), "pairwise", adjust = "tukey"))
  
  print(summary(obj.lmer6))
  print(anova(obj.lmer6, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer6, specs=c("NovelDay","Environ")), "pairwise", adjust = "tukey"))
  
  print(summary(obj.lmer7))
  print(anova(obj.lmer7, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer7, specs=c("NovelDay","Environ")), "pairwise", adjust = "tukey"))
  
  print(summary(obj.lmer8))
  print(anova(obj.lmer8, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer8, specs=c("NovelDay","Environ")), "pairwise", adjust = "tukey"))
  
  sink()
}
