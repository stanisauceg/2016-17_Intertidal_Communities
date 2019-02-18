library(tidyverse)
library(vegan)

### Ordination ###

### nMDS #########

names(data)
community <- data[,6:16]
id.tags <- data[,c(1:4,17:18)]

# which observations are for which times, given sequence 0 - 1 - 4 - 0 - 1 - 4 - etc.
t0 <- seq(from=1,by=3,length=108)
t1 <- seq(from=2,by=3,length=108)
t4 <- seq(from=3,by=3,length=108)

# removal treatment for each observation
treat <- id.tags$removal[t0]


### run nMDS 
set.seed(24)
community.mds<-metaMDS(community,distance="bray",trymax=50,k=3)

community.mds

## check stress plot
stressplot(community.mds)

## create data frame to plot nMDS with ggplot2
NMDS <- tibble(MDS1 = community.mds$points[,1], 
               MDS2 = community.mds$points[,2], 
               MDS3 = community.mds$points[,3],
               Plot = data$plot, 
               Site = data$site, 
               TimeStep = data$TimeStep,
               Removal = data$removal,
               Pred = data$pred)
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
mds_cols <- cbind(MDS1 = NMDS$MDS1, MDS2 = NMDS$MDS2, MDS3 = NMDS$MDS3)

centroids <- aggregate(mds_cols ~ Removal * Pred * TimeStep, NMDS, mean) # full factorial, 3 factors
centroids.r <- aggregate(mds_cols ~ Removal * TimeStep, NMDS, mean) # no predation

rm(mds_cols)

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


##### Analysis of Similarities (ANOSIM) ###################

my_anosims <- function(timestep) {
  a <- anosim(data[timestep,6:16], grouping=data$removal[timestep])
  print(summary(a))
  plot(a)
}

list(t0, t1, t4) %>% walk(my_anosims)


##### PERMANOVA ##########################

# full model
set.seed(24)
adonis.full <- adonis(vegdist(data[,c(5:16)], method="bray")
                    ~ TimeStep*pred*removal, strata=data$site, data=data, permutations=999)
adonis.full

adonis.2 <- adonis(vegdist(data[,c(5:16)], method="bray")
                      ~ TimeStep*pred*removal -TimeStep:pred:removal , strata=data$site, data=data, permutations=999)
adonis.2


# split by timestep
adonis_time <- function(timestep) {
  a <- adonis(vegdist(data[timestep, c(6:16)], method="bray")
              ~ pred*removal, strata=data$site[timestep], data=data[timestep,], permutations=9999)
  return(a)
}

adonis.t0 <- adonis_time(t0)
adonis.t1 <- adonis_time(t1)
adonis.t4 <- adonis_time(t4)

adonis.t0
adonis.t1
adonis.t4
