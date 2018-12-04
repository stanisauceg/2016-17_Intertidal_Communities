getwd()
setwd("C:/Users/Misha/Dropbox/research/2016-17 intertidal/")

search() # check for attached objects, detach() anything unnecessary
rm(list=ls())

library(data.table)
data<-fread("cover-wt.csv",header=TRUE)
data[,1]<-NULL
colnames(data)[8]<-"fraction"

data$target<-as.factor(data$target)
library(car)

lm1 <- lm(iA ~ target + final + iB + iM, data = data)
summary(lm1)

anova.lm1<-anova(lm1)
aov.lm1<-aov(lm1)
anova.lm1
aov.lm1
summary(aov.lm1)
summary(anova.lm1)

avPlots(lm1)
plot(lm1)