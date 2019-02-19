### create data frame showing cover of each target species over time
# melt:
mdat <- data.table::melt(data, id.vars= c("site","plot","TimeStep","date","pred","removal"),
             measure.vars = c("Fucus","Mytilus","barnacle"), 
             variable.name="species", value.name="pct.cover")
head(mdat)
summary(mdat)
str(mdat)

mdat$sp_t <- paste(mdat$species, mdat$TimeStep, sep="-")
mdat$sp_t <- factor(mdat$sp_t)

species_names <- c("barnacle" = "Barnacle", 
                  "Fucus" = "Fucus", 
                  "Mytilus" = "Mytilus", 
                  "Control" = "Control")

dodge <- position_dodge(width = 0.1)

ggplot(mdat, aes(x = TimeStep, y = pct.cover, color = removal, pch = pred)) +
  stat_summary(fun.y = mean, geom = "point", size = 2, position = dodge) +
  stat_summary(fun.y = mean, geom = "line", aes(group = interaction(removal, pred)), 
               linetype = 2, position = dodge) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1, position = dodge) +
  facet_wrap(~species, nrow = 3, labeller = as_labeller(species_names)) +
  scale_x_discrete(name = "Time Step", labels = c("Initial", "Treatment", "Recovery")) +
  ylab("Percent Cover") +
  scale_color_manual(name="Removal", labels = c("Fucus", "Barnacle", "Mytilus", "Control"),
                     values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF")) +
  theme_bw(base_size = 18)
