library(car)
library(tidyverse)

path.read.csv <- function(...){
  data_path <- file.path(...)
  df <- read_csv(data_path)
  return(df)
}

scrapings <- path.read.csv("data", "cover-wt.csv")

colnames(scrapings)[9]<-"fraction"

scrapings$target <- as.factor(scrapings$target)

lm1 <- lm(iA ~ target + final + iB + iM, data = scrapings)
summary(lm1)

anova.lm1 <- anova(lm1)
anova.lm1

# aov.lm1 <- aov(lm1)
# aov.lm1
# summary(aov.lm1)

avPlots(lm1)
plot(lm1)