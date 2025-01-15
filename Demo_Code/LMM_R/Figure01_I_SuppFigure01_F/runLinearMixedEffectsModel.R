################################################################################
################################################################################
############### RAMP-DOWN ACTIVITY PER POSITION BIN#############################
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
dt = "resFR_allBins_timebin"

tblFilename = paste0(foldN, "/TableData_", dt, ".txt")
txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create .txt file to save LMM outputs

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
mydata$StimCond = as.factor(mydata$StimCond)

sink(txtFilename)

for (ct in c(1:3)) {
  dat <- mydata[is.element(mydata$AnimalID, wtMice) &
                  mydata$CellType == ct &
                  mydata$Environ == 1, ]
  obj.lmer = lmerTest::lmer(DataY ~ PosBin + (1|AnimalID/CellID), data = dat)

  print(paste0(cellTypes[ct], " Familiar AllDaysCombined"))
  print(summary(obj.lmer))
  print(anova(obj.lmer))
  print(contrast(emmeans(obj.lmer, specs=c("PosBin")), "trt.vs.ctrl", adjust = "sidak"))
}
sink()

rm()
dt = "resFR_allBins_position"

tblFilename = paste0(foldN, "/TableData_", dt, ".txt")
txtFilename = paste0(foldN, "/stats_lme4_", dt, ".txt") #create .txt file to save LMM outputs

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
mydata$StimCond = as.factor(mydata$StimCond)

sink(txtFilename)

for (ct in c(1:3)) {
  dat <- mydata[is.element(mydata$AnimalID, wtMice) &
                  mydata$CellType == ct &
                  mydata$Environ == 1, ]
  obj.lmer = lmerTest::lmer(DataY ~ PosBin + (1|AnimalID/CellID), data = dat)

  print(paste0(cellTypes[ct], " Familiar AllDaysCombined"))
  print(summary(obj.lmer))
  print(anova(obj.lmer))
  print(contrast(emmeans(obj.lmer, specs=c("PosBin")), "trt.vs.ctrl", adjust = "sidak"))
}
sink()
