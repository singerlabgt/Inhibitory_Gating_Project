# 09102024: XZ performs session wise SWR analysis by averaging SWR duration, power, rate by session. 
### load the libraries
library(lme4)
library(emmeans)
library(sjPlot)
library(report)
library(ggplot2)

wtMice <- c(11,18,21,24,45,46,47)
pvGSMice <- c(57,61,62,63,65,4)
pvFamStimMice <- c(48,50,52,53,54,57,61,62,63,65,4)
cellTypes <-c("PV","Pyr") #corresponds to 1-3

# ################################################################################
# ################################################################################
# ############### RAMP-DOWN ACTIVITY PER POSITION BIN#############################
# 
foldN = "."

## lme4 on reward-centered firing rates
# pbkrtest.limit = 12769
# lmerTest.limit = 12769
tblFilename = paste0(foldN, "/TableData_resFR_allBins_timebin.txt")
txtFilename = paste0(foldN, "/stats_lme4_resFR_allBins_timebin.txt") #create .txt file to save LMM outputs

mydata <- read.csv(tblFilename, head=T)
mydata <- mydata[!is.na(mydata$DataY),]

#factorize variables to be treated as categories/levels instead of numerical value
mydata$AnimalID = as.factor(mydata$AnimalID)
mydata$RecDate = as.factor(mydata$RecDate)
mydata$NovelDay = as.factor(mydata$NovelDay)
mydata$PosBin = as.factor(mydata$PosBin)
mydata$Environ = as.factor(mydata$Environ)
mydata$CellID = as.factor(mydata$CellID)
mydata$CellType = as.factor(mydata$CellType)


sink(txtFilename)
for (iNovelDay in unique(mydata$NovelDay)){
  for (iCT in c(1, 2)){
    dat <- mydata[mydata$CellType == iCT & mydata$NovelDay == iNovelDay,]
    obj.lmer = lmerTest::lmer(DataY ~ PosBin + (1|AnimalID/CellID), data = dat)
    
    print(paste0(cellTypes[iCT], " NOV Day ", iNovelDay))
    print(summary(obj.lmer))
    print(anova(obj.lmer))
    print(contrast(emmeans(obj.lmer, specs=c("PosBin")), "trt.vs.ctrl", adjust = "sidak"))
  }
}

sink()
rm()