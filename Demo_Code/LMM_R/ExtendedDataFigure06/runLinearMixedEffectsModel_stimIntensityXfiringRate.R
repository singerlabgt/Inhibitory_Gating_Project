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
# 
foldN = "."
celltypes = c("PV", "pyr")
analysistypes = c("HighStim", "lowStim", "Nostim", "WT")
# process no-SWR firing rates. interested in main positional effect
pbkrtest.limit = 21320
lmerTest.limit = 21320
for (ct in 1:length(celltypes)){
  celltype = celltypes[ct]
  for (i in 1:length(analysistypes)){
    analysistype = analysistypes[i]
    tblFilename = paste0(foldN, "/TableData_resFR_", celltype, "_", analysistype, ".txt")
    txtFilename = paste0(foldN, "/stats_lme4_resFR_", celltype, "_", analysistype, ".txt") #create .txt file to save LMM outputs
    
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
    
    
    dat <- mydata
    obj.lmer = lmerTest::lmer(DataY ~ PosBin + (1|AnimalID/CellID), data = dat)
    
    print(paste0(celltype, " NOV AllDaysCombined"))
    print(summary(obj.lmer))
    print(anova(obj.lmer))
    print(contrast(emmeans(obj.lmer, specs=c("PosBin")), "trt.vs.ctrl", adjust = "sidak"))
    
    sink()
    rm()
  }
  
}

