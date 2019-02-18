library(tidyverse)
library(vegan)

# filter out four plots that were scraped before a pre-treatment photo could be taken
data_del <- data %>% 
  filter(!plot %in% c("L2-19", "L2-23", "L2-24", "L2-26") | TimeStep != "0")

summary(data_del)
names(data_del)

# split into community data & environment data (id_tags)
community <- data_del[,6:16]
id_tags <- data_del[,c(1:4,17:18)]

t0 <- c(which(data_del$TimeStep==0))
t1 <- c(which(data_del$TimeStep==1))
t4 <- c(which(data_del$TimeStep==4))

set.seed(24)
community_mds<-metaMDS(community,distance="bray",trymax=50,k=3)
### check stress plot
stressplot(community_mds)

NMDS <- data.frame(MDS1 = community_mds$points[,1], 
                   MDS2 = community_mds$points[,2], 
                   MDS3 = community_mds$points[,3], 
                   Plot = data_del$plot, 
                   Site = data_del$site, 
                   TimeStep = data_del$TimeStep,
                   Removal = data_del$removal, 
                   Pred = data_del$pred)


adonis_del_time <- function(timestep) {
  a <- adonis(vegdist(data_del[timestep, c(6:16)], method="bray")
         ~ pred*removal, strata=data_del$site[timestep], data=data_del[timestep,], permutations=9999)
  return(a)
}

adonis_del_time(t0)
adonis_del_time(t1)
adonis_del_time(t4)

# yields (qualitatively) same PERMANOVA results as unrestricted data set