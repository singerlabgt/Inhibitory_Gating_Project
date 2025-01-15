### load the libraries
library(lme4)
library(emmeans)
library(sjPlot)
library(report)
library(ggplot2)

wtMice <- c(11,18,21,24,45,46,47)
pvGSMice <- c(57,61,62,63,65,4)
pvFamStimMice <- c(48,50,52,53,54,57,61,62,63,65,4)
cellTypes <-c("NarrowInt","WideInt","Pyr") #corresponds to 1-3


foldN = "."

################################################################################
################################################################################
######################## stim effects on FR ####################################
dataStrings = c("rawFRxIntensity","normFRxIntensity")
cellTypes <-c("PV", "NarrowInt","WideInt","Pyr") #corresponds to 1-4

for (dt in c(1,2)) {
  tblFilename = paste0(foldN, "/TableData_", dataStrings[dt], ".txt")
  txtFilename = paste0(foldN, "/stats_lme4_", dataStrings[dt], ".txt") #create .txt file to save LMM outputs 
  
  mydata <- read.csv(tblFilename, head=T)
  mydata <- mydata[!is.na(mydata$DataY), ]
  
  #factorize variables to be treated as categories/levels instead of numerical value
  mydata$AnimalID = as.factor(mydata$AnimalID)
  mydata$NovelDay = as.factor(mydata$NovelDay)
  mydata$Environ = as.factor(mydata$Environ)
  mydata$CellID = as.factor(mydata$CellID)
  mydata$CellType = as.factor(mydata$CellType)
  mydata$StimCond = as.factor(mydata$StimCond)
  mydata$StimIntensity = as.factor(mydata$StimIntensity)
  
  sink(txtFilename)
  for (ct in c(1:4)) {
    pv_fam_ct <- mydata[is.element(mydata$AnimalID, pvFamStimMice) & mydata$Environ == 1 & is.element(mydata$StimIntensity, c(1,2,3)) & mydata$CellType == ct, ] 
    obj.lmer = lmerTest::lmer(DataY ~ StimIntensity + (1|AnimalID), data = pv_fam_ct)
    print(paste0(cellTypes[ct], " - Fam"))
    print(summary(obj.lmer))
    print(anova(obj.lmer))
    print(contrast(emmeans(obj.lmer, specs=c("StimIntensity")), "pairwise", adjust = "tukey"))
    
    
    pv_nov_ct <- mydata[is.element(mydata$AnimalID, pvGSMice) & mydata$Environ == 2 & is.element(mydata$StimIntensity, c(1,2,3)) & is.element(mydata$StimCond, c(2,3)) & mydata$CellType == ct, ] 
    print(paste0(cellTypes[ct], " - Nov"))
    obj.lmer = lmerTest::lmer(DataY ~ StimCond*StimIntensity + (1|AnimalID/CellID), data = pv_nov_ct)
    print(summary(obj.lmer))
    print(anova(obj.lmer))
    print(contrast(emmeans(obj.lmer, specs=c("StimIntensity","StimCond")), "pairwise", adjust = "tukey"))
  }
  sink()
}

