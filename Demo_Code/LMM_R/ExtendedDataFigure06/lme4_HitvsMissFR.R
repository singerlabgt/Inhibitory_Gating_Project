
### load the libraries
rm(list = setdiff(ls(), lsf.str()))

library(nlme)
library(lme4)
library(lmerTest)
library(emmeans)

foldN = "."

### 

tblFilename = paste0(foldN, "/data_for_lmmR_fam.txt")
txtFilename = paste0(foldN, "/stats_lme4_fam.txt") #create .txt file to save LMM outputs 

wtMice <- c(11,18,21,24,45,46,47)
pvGSMice <- c(57,61,62,63,65,4)
pvFamStimMice <- c(48,50,52,53,54,57,61,62,63,65,4)
dataStrings <- c("spatialInfo","ratemapCorr")
cellTypes <-c("NarrowInt","Pyr") #corresponds to 1-2

famdata = read.csv(tblFilename, head=T)
famdata$AnimalID = as.factor(famdata$AnimalID)
famdata$NovelDay = as.factor(famdata$NovelDay)
famdata$isCorrect = as.factor(famdata$isCorrect)
famdata$CellID = as.factor(famdata$CellID)
famdata$CellType = as.factor(famdata$CellType)

nsInt <- famdata[is.element(famdata$AnimalID, wtMice) & famdata$CellType==1 , ] 
Pyr <- famdata[is.element(famdata$AnimalID, wtMice) & famdata$CellType==2 , ] 

obj.lmer1 <- lmerTest::lmer(resFR ~ isCorrect + (1|AnimalID/CellID), data=nsInt)
obj.lmer3 <- lmerTest::lmer(resFR ~ isCorrect + (1|AnimalID/CellID), data=Pyr)

sink(txtFilename)
print(summary(obj.lmer1))
print(anova(obj.lmer1, ddf="Kenward-Roger"))
print(contrast(emmeans(obj.lmer1, specs=c("isCorrect")), "pairwise", adjust = "tukey"))

print(summary(obj.lmer3))
print(anova(obj.lmer3, ddf="Kenward-Roger"))
print(contrast(emmeans(obj.lmer3, specs=c("isCorrect")), "pairwise", adjust = "tukey"))
sink()

###

tblFilename = paste0(foldN, "/data_for_lmmR_nov.txt")
txtFilename = paste0(foldN, "/stats_lme4_nov.txt") #create .txt file to save LMM outputs 


novdata = read.csv(tblFilename, head=T)
novdata$AnimalID = as.factor(novdata$AnimalID)
novdata$NovelDay = as.factor(novdata$NovelDay)
novdata$isCorrect = as.factor(novdata$isCorrect)
novdata$CellID = as.factor(novdata$CellID)
novdata$CellType = as.factor(novdata$CellType)

nsInt <- novdata[is.element(novdata$AnimalID, wtMice) & novdata$CellType==1 , ] 
Pyr <- novdata[is.element(novdata$AnimalID, wtMice) & novdata$CellType==2 , ] 

obj.lmer1 <- lmerTest::lmer(resFR ~ isCorrect + (1|AnimalID/CellID), data=nsInt)
obj.lmer3 <- lmerTest::lmer(resFR ~ isCorrect + (1|AnimalID/CellID), data=Pyr)

sink(txtFilename)
print(summary(obj.lmer1))
print(anova(obj.lmer1, ddf="Kenward-Roger"))
print(contrast(emmeans(obj.lmer1, specs=c("isCorrect")), "pairwise", adjust = "tukey"))

print(summary(obj.lmer3))
print(anova(obj.lmer3, ddf="Kenward-Roger"))
print(contrast(emmeans(obj.lmer3, specs=c("isCorrect")), "pairwise", adjust = "tukey"))
sink()



