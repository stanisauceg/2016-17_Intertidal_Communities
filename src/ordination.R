library(tidyverse)
library(vegan)

### Ordination ###

### nMDS #########

names(data)
community <- data[,6:16]
id.tags <- data[,c(1:4,17:18)]

# which observations are for which times, given sequence 0 - 1 - 4 - 0 - 1 - 4 - etc.
t0<-seq(from=1,by=3,length=108)
t1<-seq(from=2,by=3,length=108)
t4<-seq(from=3,by=3,length=108)

# removal treatment for each observation
treat<-id.tags$removal[t0]


### run nMDS 
set.seed(24)
community.mds<-metaMDS(community,distance="bray",trymax=50,k=3)

## check stress plot
stressplot(community.mds)


# create data frame to plot nMDS with ggplot2
MDS1<-community.mds$points[,1]
MDS2<-community.mds$points[,2]
MDS3<-community.mds$points[,3]

NMDS <- tibble(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3, 
                   Plot = data$plot, Site = data$site, TimeStep = data$TimeStep,
                   Removal = data$removal, Pred = data$pred)
summary(NMDS)
str(NMDS)

# Rename removal treatments to be more informative. 
# Currently a factor, so un-factor it first.
NMDS$Removal <- as.character(NMDS$Removal)

# rename treatments
NMDS$Removal[which(NMDS$Removal == "A")] <- "Fucus"
NMDS$Removal[which(NMDS$Removal == "B")] <- "Barnacle"
NMDS$Removal[which(NMDS$Removal == "M")] <- "Mussel"
NMDS$Removal[which(NMDS$Removal == "none")] <- "No Removal"

# re-factor
NMDS$Removal <- factor(NMDS$Removal, levels = c("Fucus", "Barnacle", "Mussel", "No Removal"))

# get the centroids for each grouping
mds_cols <- cbind(MDS1, MDS2, MDS3)

centroids <- aggregate(mds_cols ~ Removal * Pred * TimeStep, NMDS, mean) # full factorial, 3 factors
centroids.r <- aggregate(mds_cols ~ Removal * TimeStep, NMDS, mean) # no predation

# species nMDS scores
NMDS_species <- data.frame(scores(community.mds, display="species"))
rownames(NMDS_species)
rownames(NMDS_species)[5]<-"encrusting"
rownames(NMDS_species)[7]<-"debris"


#### NMDS plots ###########

# show ggplot default colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

show_colors <- function(n) {
  cols = gg_color_hue(n)
  plot(1:n, pch = 16, cex = 2, col = cols)
}

show_colors(4)

# show ggplot default line types
show_lines <- function(n) {
  lines <- data.frame(x = 1:2, y = rep(n:1, each = 2), grp = factor(rep(1:n, each = 2)))
  
  ggplot(data = lines, aes(x = x, y = y, linetype = grp)) +
    geom_line() +
    geom_text(aes(x = 0.95, label = grp)) +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = "none")
}

show_lines(13)

### Site groupings ----

site_plot <- function(timestep, plot_title) {
  g <- ggplot(NMDS[timestep,], aes(x=MDS1, y=MDS2, col=Site)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    labs(title = plot_title)+
    coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))+
    geom_text(inherit.aes=FALSE, data=NMDS_species, label=rownames(NMDS_species), 
              aes(x=NMDS1, y=NMDS2), color="red")
  print(g)
}

site_plot(t0, "Initial (pre-treatment)")

### Removal treatment groupings ----

removal_plot <- function(timestep, plot_title){
  g <- ggplot(NMDS[timestep,], aes(x = MDS1, y = MDS2, col = Removal)) +
    geom_point() +
    stat_ellipse() +
    scale_color_manual(values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF")) +
    theme_bw() +
    labs(title = plot_title) +
    coord_cartesian(xlim = c(-1,1), ylim = c(-1,1)) +
    geom_text(inherit.aes = FALSE, data = NMDS_species[c(1,2,4),], label = rownames(NMDS_species)[c(1,2,4)], 
              aes(x = NMDS1, y = NMDS2), color = "red")
  print(g)
}

removal_plot(t0, "Initial (pre-treatment)")


### Predator groupings ####

pred_plot <- function(timestep, plot_title){
  g <- ggplot(NMDS[timestep,], aes(x=MDS1, y=MDS2, col=Pred)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    labs(title = plot_title)
  print(g)
}

pred_plot(t0, "Initial (pre-treatment)")

### Pred * Removal groupings ####

pred_by_removal_plot <- function(timestep, plot_title) {
  
  num_step <- NMDS[timestep,]$TimeStep %>% as.character() %>% as.numeric() %>% mean()
  
  g <- ggplot(NMDS[timestep,], aes(x=MDS1, y=MDS2, col=Removal, pch=Pred)) +
    scale_color_manual(values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF")) +
    scale_linetype_manual(values = c(1,3,2)) +
    stat_ellipse(aes(linetype = Pred)) +
    theme_bw() +
    labs(title = plot_title) +
    geom_point(data = centroids[which(centroids$TimeStep == num_step), ], size = 2) +
    coord_cartesian(xlim = c(-1,1.1), ylim = c(-1.1,1.5)) +
    geom_text(inherit.aes = FALSE, data = NMDS_species, label = rownames(NMDS_species), 
              aes(x = NMDS1, y = NMDS2), color = "red")
  print(g)
}

pred_by_removal_plot(t0, "Initial (pre-treatment)")

### multi-plot function to plot t0, t1, and t4 ####

plot_3_timesteps <- function(plot_fxn) {
  plot_titles <- list("Initial (pre-treatment)", "Treatment", "Recovery")
  time_steps <- list(t0, t1, t4)
  
  walk2(.f = plot_fxn, .x = time_steps, .y = plot_titles)
}

plot_3_timesteps(site_plot)
plot_3_timesteps(pred_plot)
plot_3_timesteps(removal_plot)
# plot_3_timesteps(pred_by_removal_plot)


community.mds

### Analysis of Similarities (ANOSIM) ####

my_anosims <- function(timestep) {
  a <- anosim(data[timestep,6:16], grouping=data$removal[timestep])
  print(summary(a))
  plot(a)
}

list(t0, t1, t4) %>% walk(my_anosims)


### PERMANOVA

# full model
set.seed(24)
adonis.full <- adonis(vegdist(data[,c(5:16)], method="bray")
                    ~ TimeStep*pred*removal, strata=data$site, data=data, permutations=999)

# split by timestep
t<-t0 # choose from c(t0, t1, t4)

adonis.t0<-adonis(vegdist(data[t,c(6:16)],method="bray")
                  ~pred*removal, strata=data$site[t], data=data[t], permutations=9999)
adonis.t0m<-adonis(vegdist(data[t,c(6:8,10:16)],method="bray")
                   ~pred*removal, strata=data$site[t], data=data[t], permutations=9999)
adonis.t0b<-adonis(vegdist(data[t,c(6,8:16)],method="bray")
                   ~pred*removal, strata=data$site[t], data=data[t], permutations=9999)
adonis.t0f<-adonis(vegdist(data[t,c(7:16)],method="bray")
                   ~pred*removal, strata=data$site[t], data=data[t], permutations=9999)

adonis.full
adonis.t0
adonis.t0$aov.tab
adonis.t0$aov.tab[6]
adonis.t0$aov.tab[6][[1]]
adonis.t0$aov.tab[6][[1]][[3]] # interaction p-val
adonis.t0$aov.tab[5][[1]][[3]] # interaction r2
adonis.t1
adonis.t4

adonis.t4a
adonis.t4b
adonis.t4m
adonis.t4bm

adonis.t0s
adonis.t1s
adonis.t4s

adonis2.t0<-adonis2(vegdist(data[t,c(6:16)],method="bray")
                    ~pred*removal, strata=data$site[t], data=data[t], permutations=9999)

adonis2.t1<-adonis2(vegdist(data[t,c(6:16)],method="bray")
                    ~pred*removal, strata=data$site[t], data=data[t], permutations=9999)

adonis2.t4<-adonis2(vegdist(data[t,c(6:16)],method="bray")
                    ~pred*removal, strata=data$site[t], data=data[t], permutations=9999)


adonis2.t0
adonis2.t1
adonis2.t4

dbrda.t0 <- dbrda(vegdist(data[t,c(6:16)],method="bray")
                  ~pred*removal, strata=data$site[t], data=data[t])
anova(dbrda.t0, by = "margin")


rmv<-c("none","A","B","M")
prd<-c("open","roof","cage")
combo<-expand.grid(rmv, prd)
colnames(combo)[1]<-"rmv"
colnames(combo)[2]<-"prd"

data.0<-data[t0]

data.0$permanova.rmplot<-numeric(length=length(t0))
for(i in 1:length(data.0$plot)){
  adonis.rmplot<-adonis(vegdist(data.0[-i,c(6:16)],method="bray")
                        ~pred*removal, strata=data.0$site[-i], data=data.0[-i,], permutations=999)
  data.0$permanova.rmplot[i]<-adonis.rmplot$aov.tab[6][[1]][[3]]
}
data.0$permanova.rmplot
sort(data.0$permanova.rmplot)
which(data.0$permanova.rmplot==max(data.0$permanova.rmplot))
which(data.0$permanova.rmplot==max(data.0$permanova.rmplot[-c(31,50,49,98)]))
which(data.0$permanova.rmplot==max(data.0$permanova.rmplot[-c(18,50)]))
# in descending order: 50, 18, 20 & 103 (tied)
# another time: 31 & 98 (tied), 28 & 59 (tied)
# last time: 46 & 98 (tied), 59 & 63 & 67
data.0[which(removal=="none" & pred =="cage"),2:9][order(Mytilus, decreasing=TRUE)]
data.0[which(removal=="none" & pred =="cage"),2:9][order(barnacle, decreasing=FALSE)]
which(data.0$plot=="L2-14")
which(data.0$plot=="L3-22")
which(data.0$plot=="L2-18")
which(data.0$plot=="L3-2")
which(data.0$plot=="L3-18")


i<-c(76,75,1,51,69,18) # six plots from none/roof
i<-c(42, 87, 46, 84, 82) # 5 plots from none/cage

adonis(vegdist(data.0[-i,c(6:16)],method="bray")
       ~pred*removal, strata=data.0$site[-i], data=data.0[-i,], permutations=9999)
mean(data.0$Fucus[i])
mean(data.0$Fucus[-i])
mean(data.0$Mytilus[i])
mean(data.0$Mytilus[-i])
mean(data.0$barnacle[i])
mean(data.0$barnacle[-i])
sort(data.0$plot[i])

data.4<-data[t4]
i<-which((data.4$removal=="M" | data.4$removal=="B") & data.4$pred=="open") # six plots from none/roof
adonis(vegdist(data.4[,c(6:16)],method="bray")
       ~pred*removal, strata=data.4$site, data=data.4, permutations=9999)
adonis(vegdist(data.4[-i,c(6:16)],method="bray")
       ~pred*removal, strata=data.4$site[-i], data=data.4[-i,], permutations=9999)


impt<-which(data.0$permanova.rmplot>0.013)
data.0$removal[impt]
data.0$pred[impt]
data.0$plot[impt]

permanova.rmcombo<-numeric(length=length(combo[,1]))
r2.rmcombo<-numeric(length=length(combo[,1]))
for (i in 1:length(permanova.rmcombo)) {
  sp<-combo[i,1]
  pred<-combo[i,2]
  data.c<-data.0[-which(data.0$removal==sp & data.0$pred==pred),]
  adonis.rmcombo<-adonis(vegdist(data.c[,c(6:16)],method="bray")
                         ~pred*removal, strata=data.c$site, data=data.c, permutations=999)
  permanova.rmcombo[i] <- adonis.rmcombo$aov.tab[6][[1]][[3]]
  r2.rmcombo[i] <- adonis.rmcombo$aov.tab[5][[1]][[3]]
}
permanova.rmcombo
r2.rmcombo
impt<-which(permanova.rmcombo>0.05)
which(permanova.rmcombo > 0.05)
r2.rmcombo[which(permanova.rmcombo > 0.05)]

combo[c(1,5,9),]
which(data.0$permanova.rmplot==max(data.0$permanova.rmplot))
data.0$plot[50]
data.0$removal[50]
data.0$pred[50]
data.0$plot[which(data.0$removal=="none" & data.0$pred=="cage")]


adonis.rmplot<-adonis(vegdist(data.0[-i,c(6:16)],method="bray")
                      ~pred*removal, strata=data.0$site[-i], data=data.0[-i,], permutations=999)

adonis(vegdist(data.0[,c(6:16)],method="bray")~pred*removal, strata=data.0$site,data=data.0,permutations=999)

data.no<-data.0[-which(data.0$removal=="none" & data.0$pred=="open"),]
adonis(vegdist(data.no[,c(6:16)],method="bray")
       ~pred*removal, strata=data.no$site, data=data.no, permutations=999) # still sig
data.nr<-data.0[-which(data.0$removal=="none" & data.0$pred=="roof"),]
adonis(vegdist(data.nr[,c(6:16)],method="bray")
       ~pred*removal, strata=data.nr$site, data=data.nr, permutations=999) # sig/ marginally sig
data.nc<-data.0[-which(data.0$removal=="none" & data.0$pred=="cage"),]
adonis(vegdist(data.nc[,c(6:16)],method="bray")
       ~pred*removal, strata=data.nc$site, data=data.nc, permutations=999) # still sig

adonis(vegdist(data.4[,c(6:16)],method="bray")
       ~pred*removal, strata=data.4$site, data=data.4, permutations=9999)
i<-which((data.4$removal=="M" | data.4$removal == "B") & data.4$pred=="open") # 18 plots from (-M|-B) & open
adonis(vegdist(data.4[-i,c(6:16)],method="bray")
       ~pred*removal, strata=data.4$site[-i], data=data.4[-i,], permutations=9999) # test w/o the above plots
adonis(vegdist(data.4[,c(6:8,10:16)],method="bray")
       ~pred*removal, strata=data.4$site, data=data.4[,-9], permutations=9999) # test w/o mussels in community



adonis(vegdist(data.mo[,c(6:16)],method="bray")
       ~pred*removal, strata=data.mo$site, data=data.mo, permutations=999) # still sig


