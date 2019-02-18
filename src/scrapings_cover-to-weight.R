## estimating initial algal cover for "L2-19", "L2-23", "L2-24", and barnacle cover for "L2-26"

# L2-26 initial barnacle cover estimatable from first set of post-scrape photos (scars/shading on rock)
#   # --> best estimate is 31% barnacle cover scraped bare, 33% total initial barnacle (but canopy layer unknown!) 
#   #     high-barnacle-cover sample all fit into cup, so assuming < 1% mussel/algae removed!!

# remaining estimates written down memory af the end of the day when the plots were scraped:
# L2-19 initially ~35% Fucus cover
# L2-23 initial Fucus cover unknown - but low in surroundings, estimate 10%
# L2-24 ~90% Fucus cover

# after adding algal cover, reduce cover of remaining taxa proportionally so that total cover stays 100%

library(tidyverse)

path.read.csv <- function(...){
  data_path <- file.path(...)
  df <- read_csv(data_path)
  return(df)
}

scrapings <- path.read.csv("data", "cover-wt-2.csv")

scrapings <- scrapings %>%
  mutate(removal = as.factor(removal)) %>%
  rename(iA = Fucus,
         iB = barnacle,
         iM = Mytilus)

fucus_removals <- scrapings %>%
  filter(removal == "A")

# version 1: assume linear relationship btw % cover of a focal taxon and mass removed
# so, final_weight (i.e., net wt removed) ~ initial cover
# and initial cover ~ final_weight

lm1 <- lm(iA ~ final_weight, data = fucus_removals)
summary(lm1)

# not great: R-squared = 0.19

predict_A <- function(final_weight) {
  lm1[[1]][[1]] + final_weight * lm1[[1]][[2]]
}

weights <- c(5.96, 10.62, 36.21)
predict_A(weights)
# compare to 35, 10, and 85-90 % cover

