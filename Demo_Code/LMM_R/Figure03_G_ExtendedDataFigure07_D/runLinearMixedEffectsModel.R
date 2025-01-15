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

# ################################################################################
# ################################################################################
# ############### RAMP-DOWN ACTIVITY PER POSITION BIN#############################
foldN = "."

# process no-SWR firing rates. interested in main positional effect
tblFilename = paste0(foldN, "/TableData_resFR_allBins_FirstLastBlock_PV.txt")
txtFilename = paste0(foldN, "/stats_lme4_resFR_allBins_FirstLastBlock_PV.txt") #create .txt file to save LMM outputs

mydata <- read.csv(tblFilename, head=T)
mydata <- mydata[!is.na(mydata$DataY),]

#factorize variables to be treated as categories/levels instead of numerical value
mydata$TrialBlock = as.factor(mydata$TrialBlock)
mydata$PosBin = as.factor(mydata$PosBin)

sink(txtFilename)

for (iBlock in c(1, 2)){
  dat <- mydata[mydata$TrialBlock == iBlock,]
  obj.lm = lm(DataY ~ PosBin, data = dat)
  print(paste0("PV cells, Nov Day 1, Trial block ", iBlock))
  print(summary(obj.lm))
  print(anova(obj.lm))
  print(anova(obj.lmer))
  print(contrast(emmeans(obj.lm, specs=c("PosBin")), "trt.vs.ctrl", adjust = "sidak"))
}


sink()
rm()
################################################################################
# Pyramidal cell spatial info
tblFilename = paste0(foldN, "/TableData_SpatialInfo_allTrialBlocks_pyr.txt")
txtFilename = paste0(foldN, "/stats_lme4_SpatialInfo_allTrialBlocks_pyr.txt") #create .txt file to save LMM outputs

mydata <- read.csv(tblFilename, head=T)
mydata <- mydata[!is.na(mydata$DataY),]

#factorize variables to be treated as categories/levels instead of numerical value
mydata$TrialBlock = as.factor(mydata$TrialBlock)
mydata$AnimalID = as.factor(mydata$AnimalID)
mydata$CellID = as.factor(mydata$CellID)

sink(txtFilename)


dat <- mydata
obj.lmer = lmerTest::lmer(DataY ~ TrialBlock + (1|AnimalID/CellID), data = dat)

print(paste0("Pyramidal cells, Nov Day 1, all trial blocks"))
print(summary(obj.lmer))
print(anova(obj.lmer))
print(contrast(emmeans(obj.lmer, specs=c("TrialBlock")), "trt.vs.ctrl", adjust = "sidak"))

sink()
rm()
################################################################################
# Pyramidal cell ratemap correlation
tblFilename = paste0(foldN, "/TableData_RatemapCorr_allTrialBlocks_pyr.txt")
txtFilename = paste0(foldN, "/stats_lme4_RatemapCorr_allTrialBlocks_pyr.txt") #create .txt file to save LMM outputs

mydata <- read.csv(tblFilename, head=T)
mydata <- mydata[!is.na(mydata$DataY),]

#factorize variables to be treated as categories/levels instead of numerical value
mydata$TrialBlock = as.factor(mydata$TrialBlock)
mydata$AnimalID = as.factor(mydata$AnimalID)
mydata$CellID = as.factor(mydata$CellID)

sink(txtFilename)


dat <- mydata
obj.lmer = lmerTest::lmer(DataY ~ TrialBlock + (1|AnimalID/CellID), data = dat)

print(paste0("Pyramidal cells, Nov Day 1, all trial blocks"))
print(summary(obj.lmer))
print(anova(obj.lmer))
print(contrast(emmeans(obj.lmer, specs=c("TrialBlock")), "trt.vs.ctrl", adjust = "sidak"))

sink()
rm()
