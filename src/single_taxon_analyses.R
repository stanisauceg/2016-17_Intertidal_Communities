### create data frame showing cover of each target species over time
mdat <- data %>% 
  gather(key = "species", value = "pct_cover", Fucus, barnacle, Mytilus) %>%
  select(TimeStep:plot, removal:pct_cover) %>%
  mutate(sp_t = factor(paste(species, TimeStep, sep = "-")))

head(mdat)
summary(mdat)
str(mdat)

# dictionary for species names
species_names <- c("barnacle" = "Barnacle", 
                  "Fucus" = "Fucus", 
                  "Mytilus" = "Mussel", 
                  "Control" = "Control")

dodge <- position_dodge(width = 0.1)

ggplot(mdat, aes(x = TimeStep, y = pct_cover, color = removal)) +
  stat_summary(fun.y = mean, geom = "point", size = 2, position = dodge) +
  stat_summary(fun.y = mean, geom = "line", aes(group = removal), 
               linetype = 2, position = dodge) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1, position = dodge) +
  facet_wrap(~ species, nrow = 3, labeller = as_labeller(species_names)) +
  scale_x_discrete(name = "Time Step", labels = c("Initial", "Treatment", "Recovery")) +
  ylab("Percent Cover") +
  scale_color_manual(name = "Removal", labels = c("Fucus", "Barnacle", "Mussel", "Control"),
                     values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF")) +
  theme_bw(base_size = 18)


# mdat2 <- data %>% 
#   gather(key = "species", value = "pct_cover", bare) %>%
#   select(TimeStep:plot, removal:pct_cover) 
# 
# ggplot(mdat2, aes(x = TimeStep, y = pct_cover, color = removal)) +
#   stat_summary(fun.y = mean, geom = "point", size = 2, position = dodge) +
#   stat_summary(fun.y = mean, geom = "line", aes(group = removal), linetype = 2, position = dodge) +
#   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1, position = dodge) +
#   xlab("Time Step") +
#   ylab("Percent Bare Space") +
#   scale_color_manual(name="Removal", labels = c("Fucus", "Barnacle", "Mytilus", "Control"),
#                      values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF")) +
#   scale_x_discrete(labels = c("Initial", "Treatment", "Recovery")) +
#   theme_bw(base_size = 18) +
#   theme(axis.text.x = element_text(hjust = 0.5, vjust = 1))
