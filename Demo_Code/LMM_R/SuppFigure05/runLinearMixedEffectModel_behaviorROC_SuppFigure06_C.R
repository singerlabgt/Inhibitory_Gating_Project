### load the libraries

rm(list = setdiff(ls(), lsf.str()))

library(nlme)
library(lme4)
library(lmerTest)
library(emmeans)

wtMice <- c(11,18,21,24,45,46,47)
pvGSMice <- c(4,57,61,62,63,65)
rocTypes <- c("Speed","Lickrate","LickLatency")
trialTypes <- c("LowStim","HighStim","NoStim","BothStim","ShortStim")


#load MATLAB output (in table format in .txt) file in R

foldN = "." #get latest folder name 

tblFilename = paste0(foldN, "/TableData_", 'behaviorDeltaROC', ".txt")
txtFilename = paste0(foldN, "/stats_lme4_", 'behaviorDeltaROC', ".txt") #create .txt file to save LMM outputs 


sink(txtFilename) ##save output data

for (behT in rocTypes ) {
  
  mydata <- read.csv(tblFilename, head=T)
  mydata <- mydata[!is.na(mydata$DataY) & mydata$BehavT == behT, ]
  
  ## Do not forget to factor the animal IDs
  mydata$AnimalID = as.factor(mydata$AnimalID)
  mydata$RecDate = as.factor(mydata$RecDate)
  mydata$NovelDay = as.factor(mydata$NovelDay)
  mydata$Environ = as.factor(mydata$Environ)
  mydata$StimCond = as.factor(mydata$StimCond)
  
  ################## Linear Mixed-effects Model ###########################
  

  #### PVxAi32 mice - Fam environment only
  ## using all trials (both stim and nostim) of all intensities 
  pv <- mydata[is.element(mydata$AnimalID, pvGSMice) & is.element(mydata$Environ, c(1)), ]
  obj.lmer4=lmerTest::lmer(DataY ~ TrialT + (1|AnimalID), data=pv)
  print(paste0("PVxAi32 - uing all trials, fam only - ", behT))
  print(summary(obj.lmer4))
  print(anova(obj.lmer4, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer4, specs="TrialT"), "pairwise", adjust = "tukey"))

  
}

sink()

