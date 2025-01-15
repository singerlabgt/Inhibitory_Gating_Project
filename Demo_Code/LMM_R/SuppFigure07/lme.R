library(lme4)
library(emmeans)
library(sjPlot)
library(report)
# peak firing rate
df <- read.csv('tabledata_peakfr.csv', sep = ',', header = F)
df$V2 <- as.factor(df$V2)
df$V4 <- as.factor(df$V4)
df$V5 <- as.factor(df$V5)
colnames(df) <- c('PeakFR', 'AnimalID', 'Date', 'NovelDay', 'CellID', 'StimCondition')
txtFilename <- "lme4result_peakfr.txt"
StimStrings <- c('GoalStim', 'ShamStim')
sink(txtFilename)
for (iStim in c(1,2)){
  print(paste0('Percentiles for StimCondition', StimStrings[iStim]) )
  print(quantile(df[df$StimCondition == iStim, "PeakFR"], probs = c(0, 0.25, 0.5, 0.75, 1)))
}
obj.lmer = lmerTest::lmer(PeakFR ~ StimCondition + (1|AnimalID/CellID), data = df)
print(summary(obj.lmer))
print(anova(obj.lmer))
print(contrast(emmeans(obj.lmer, specs=c("StimCondition")), "pairwise", adjust = "tukey"))
sink()
rm()

# vector strength
df <- read.csv('tabledata_vecstrength.csv', sep = ',', header = F)
df$V2 <- as.factor(df$V2)
df$V4 <- as.factor(df$V4)
df$V5 <- as.factor(df$V5)
colnames(df) <- c('VecStrength', 'AnimalID', 'Date', 'NovelDay', 'CellID', 'StimCondition')
txtFilename <- "lme4result_vecstrength.txt"
sink(txtFilename)
StimStrings <- c('GoalStim', 'ShamStim')
for (iStim in c(1,2)){
  print(paste0('Percentiles for StimCondition', StimStrings[iStim]) )
  print(quantile(df[df$StimCondition == iStim,"VecStrength"], probs = c(0, 0.25, 0.5, 0.75, 1)))
}
obj.lmer = lmerTest::lmer(VecStrength ~ StimCondition + (1|AnimalID/CellID), data = df)
print(summary(obj.lmer))
print(anova(obj.lmer))
print(contrast(emmeans(obj.lmer, specs=c("StimCondition")), "pairwise", adjust = "tukey"))
sink()
rm()

# theta power
df <- read.csv('tabledata_thetapower.csv', sep = ',', header = F)
df$V2 <- as.factor(df$V2)
df$V4 <- as.factor(df$V4)
df$V5 <- as.factor(df$V5)
colnames(df) <- c('ThetaPower', 'AnimalID', 'Date', 'NovelDay', 'StimCondition')
txtFilename <- "lme4result_thetapower.txt"
sink(txtFilename)
StimStrings <- c('GoalStim', 'ShamStim')
for (iStim in c(1,2)){
  print(paste0('Percentiles for StimCondition', StimStrings[iStim]) )
  print(quantile(df[df$StimCondition == iStim,"ThetaPower"], probs = c(0, 0.25, 0.5, 0.75, 1)))
}
obj.lmer = lmerTest::lmer(ThetaPower ~ StimCondition + (1|AnimalID), data = df)
print(summary(obj.lmer))
print(anova(obj.lmer))
print(contrast(emmeans(obj.lmer, specs=c("StimCondition")), "pairwise", adjust = "tukey"))
sink()
rm()