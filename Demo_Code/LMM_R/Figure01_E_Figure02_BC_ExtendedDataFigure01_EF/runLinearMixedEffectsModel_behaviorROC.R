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
foldN = "."

################################################################################
############################## Delta ROC #######################################

tblFilename = paste0(foldN, "/TableData_", 'behaviorDeltaROC', ".txt")
txtFilename = paste0(foldN, "/stats_lme4_", 'behaviorDeltaROC', ".txt") #create .txt file to save LMM outputs 


sink(txtFilename) ##save output data

for (behT in rocTypes ) {
  
  mydata <- read.csv(tblFilename, head=T)
  colnames(mydata)[5] <- 'DataY'
  mydata <- mydata[!is.na(mydata$DataY) & mydata$BehavT == behT, ]
  
  ## Do not forget to factor the animal IDs
  mydata$AnimalID = as.factor(mydata$AnimalID)
  mydata$RecDate = as.factor(mydata$RecDate)
  mydata$NovelDay = as.factor(mydata$NovelDay)
  mydata$Environ = as.factor(mydata$Environ)
  mydata$StimCond = as.factor(mydata$StimCond)
  
  ################## Linear Mixed-effects Model ###########################
  
  
  #### WT mice 
  ## specify data baed on genotype and condition of interest
  wt <- mydata[is.element(mydata$AnimalID, wtMice) & is.element(mydata$TrialT, "AllTrials"), ]
  obj.lmer1=lmerTest::lmer(DataY ~ NovelDay*Environ + (1|AnimalID), data=wt)
  print(paste0("WT MICE - uing all trials, all environments - ", behT))
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer1, specs="NovelDay","Environ"), "pairwise", adjust = "tukey"))
  
  
  wt <- mydata[is.element(mydata$AnimalID, wtMice) & is.element(mydata$TrialT, "AllTrials") & is.element(mydata$Environ, c(1,2)), ]
  obj.lmer2=lmerTest::lmer(DataY ~ NovelDay*Environ + (1|AnimalID), data=wt)
  print(paste0("WT MICE - uing all trials, fam and nov only - ", behT))
  print(summary(obj.lmer2))
  print(anova(obj.lmer2, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer2, specs="NovelDay","Environ"), "pairwise", adjust = "tukey"))
  
  
  #### PVxAi32 mice - Novel environment only
  ## using all trials (both stim and nostim) of all intensities
  pv <- mydata[is.element(mydata$AnimalID, pvGSMice) & is.element(mydata$TrialT, "AllTrials") & is.element(mydata$Environ, c(2)), ]
  obj.lmer3=lmerTest::lmer(DataY ~ NovelDay*StimCond + (1|AnimalID), data=pv)
  print(paste0("PVxAi32 - uing all trials, nov only - ", behT))
  print(summary(obj.lmer3))
  print(anova(obj.lmer3, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer3, specs="NovelDay","StimCond"), "pairwise", adjust = "tukey"))
  
  # #### PVxAi32 mice - Fam environment only
  # ## using all trials (both stim and nostim) of all intensities
  # pv <- mydata[is.element(mydata$AnimalID, pvGSMice) & is.element(mydata$Environ, c(1)), ]
  # obj.lmer4=lmerTest::lmer(DataY ~ StimCond + (1|AnimalID), data=pv)
  # print(paste0("PVxAi32 - uing all trials, fam only - ", behT))
  # print(summary(obj.lmer4))
  # print(anova(obj.lmer4, ddf="Kenward-Roger"))
  # print(contrast(emmeans(obj.lmer4, specs="StimCond"), "pairwise", adjust = "tukey"))
  
  for (tt in trialTypes ) {
    ## pv low-stim trials only in novel 
    pv <- mydata[is.element(mydata$AnimalID, pvGSMice) & mydata$TrialT==tt & is.element(mydata$Environ, c(2)), ]
    obj.lmer5=lmerTest::lmer(DataY ~ NovelDay*StimCond + (1|AnimalID), data=pv)
    print(paste0("Nov only - ", behT, " using ", tt, "Trials"))
    print(summary(obj.lmer5))
    print(anova(obj.lmer5, ddf="Kenward-Roger"))
    print(contrast(emmeans(obj.lmer5, specs="NovelDay","StimCond"), "pairwise", adjust = "tukey"))
    
    
    # ## pv only in familiar -- exclude day info here 
    # pv <- mydata[is.element(mydata$AnimalID, pvGSMice) & is.element(mydata$Environ, c(1)) & is.element(mydata$TrialT, tt) , ]
    # obj.lmer5=lmerTest::lmer(DataY ~ TrialT + (1|AnimalID), data=pv)
    # print(paste0("Fam only - ", behT, " using Stim Trials"))
    # print(summary(obj.lmer5))
    # print(anova(obj.lmer5, ddf="Kenward-Roger"))
    # print(contrast(emmeans(obj.lmer5, specs="TrialT"), "pairwise", adjust = "tukey"))
    
  }
  
}

sink()

################################################################################
############################## Raw ROC #######################################
tblFilename = paste0(foldN, "/TableData_", 'behaviorDeltaROC', ".txt") # table data have the same name as delta rocs
txtFilename = paste0(foldN, "/stats_lme4_", 'behaviorRawROC', ".txt") #create .txt file to save LMM outputs 


sink(txtFilename) ##save output data

for (behT in rocTypes ) {
  
  mydata <- read.csv(tblFilename, head=T)
  colnames(mydata)[4] <- 'DataY'
  mydata <- mydata[!is.na(mydata$DataY) & mydata$BehavT == behT, ]
  
  ## Do not forget to factor the animal IDs
  mydata$AnimalID = as.factor(mydata$AnimalID)
  mydata$RecDate = as.factor(mydata$RecDate)
  mydata$NovelDay = as.factor(mydata$NovelDay)
  mydata$Environ = as.factor(mydata$Environ)
  mydata$StimCond = as.factor(mydata$StimCond)
  
  ################## Linear Mixed-effects Model ###########################
  
  
  #### WT mice 
  ## specify data baed on genotype and condition of interest
  wt <- mydata[is.element(mydata$AnimalID, wtMice) & is.element(mydata$TrialT, "AllTrials"), ]
  obj.lmer1=lmerTest::lmer(DataY ~ NovelDay*Environ + (1|AnimalID), data=wt)
  print(paste0("WT MICE - uing all trials, all environments - ", behT))
  print(summary(obj.lmer1))
  print(anova(obj.lmer1, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer1, specs="NovelDay","Environ"), "pairwise", adjust = "tukey"))
  
  
  wt <- mydata[is.element(mydata$AnimalID, wtMice) & is.element(mydata$TrialT, "AllTrials") & is.element(mydata$Environ, c(1,2)), ]
  obj.lmer2=lmerTest::lmer(DataY ~ NovelDay*Environ + (1|AnimalID), data=wt)
  print(paste0("WT MICE - uing all trials, fam and nov only - ", behT))
  print(summary(obj.lmer2))
  print(anova(obj.lmer2, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer2, specs="NovelDay","Environ"), "pairwise", adjust = "tukey"))
  
  
  #### PVxAi32 mice - Novel environment only
  ## using all trials (both stim and nostim) of all intensities
  pv <- mydata[is.element(mydata$AnimalID, pvGSMice) & is.element(mydata$TrialT, "AllTrials") & is.element(mydata$Environ, c(2)), ]
  obj.lmer3=lmerTest::lmer(DataY ~ NovelDay*StimCond + (1|AnimalID), data=pv)
  print(paste0("PVxAi32 - uing all trials, nov only - ", behT))
  print(summary(obj.lmer3))
  print(anova(obj.lmer3, ddf="Kenward-Roger"))
  print(contrast(emmeans(obj.lmer3, specs="NovelDay","StimCond"), "pairwise", adjust = "tukey"))
  
  # #### PVxAi32 mice - Fam environment only
  # ## using all trials (both stim and nostim) of all intensities
  # pv <- mydata[is.element(mydata$AnimalID, pvGSMice) & is.element(mydata$Environ, c(1)), ]
  # obj.lmer4=lmerTest::lmer(DataY ~ StimCond + (1|AnimalID), data=pv)
  # print(paste0("PVxAi32 - uing all trials, fam only - ", behT))
  # print(summary(obj.lmer4))
  # print(anova(obj.lmer4, ddf="Kenward-Roger"))
  # print(contrast(emmeans(obj.lmer4, specs="StimCond"), "pairwise", adjust = "tukey"))
  
  for (tt in trialTypes ) {
    ## pv low-stim trials only in novel 
    pv <- mydata[is.element(mydata$AnimalID, pvGSMice) & mydata$TrialT==tt & is.element(mydata$Environ, c(2)), ]
    obj.lmer5=lmerTest::lmer(DataY ~ NovelDay*StimCond + (1|AnimalID), data=pv)
    print(paste0("Nov only - ", behT, " using ", tt, "Trials"))
    print(summary(obj.lmer5))
    print(anova(obj.lmer5, ddf="Kenward-Roger"))
    print(contrast(emmeans(obj.lmer5, specs="NovelDay","StimCond"), "pairwise", adjust = "tukey"))
    
    
    # ## pv only in familiar -- exclude day info here 
    # pv <- mydata[is.element(mydata$AnimalID, pvGSMice) & is.element(mydata$Environ, c(1)) & is.element(mydata$TrialT, tt) , ]
    # obj.lmer5=lmerTest::lmer(DataY ~ TrialT + (1|AnimalID), data=pv)
    # print(paste0("Fam only - ", behT, " using Stim Trials"))
    # print(summary(obj.lmer5))
    # print(anova(obj.lmer5, ddf="Kenward-Roger"))
    # print(contrast(emmeans(obj.lmer5, specs="TrialT"), "pairwise", adjust = "tukey"))
    
  }
  
}

sink()