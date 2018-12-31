### libraries ####

library(data.table)
library(binom)
library(vegan)
library(lme4)
library(pbkrtest)
library(tidyverse)

### load data, remove extra rows ####

data.path <- file.path("data", "counts.csv")
data <- read_csv(data.path)

dim(data)
# data should not have 2578 obs.
str(data)
head(data)
tail(data)
data[320:335,]

# remove unsorted chaff at bottom of document - this was leftover from on-the-spot calculations
data <- data[0:330,]
# still some extra rows, will deal with them soon

## load treatments
treatments.path <- file.path("data", "treatments.csv")
treatments <- read_csv(treatments.path)
dim(treatments)

head(treatments)
# retain Site & ID as key columns, diversity and pred treatment columns; only first 108 rows have useful data
treatments <- treatments[1:108,c(1,3:5)]
treatments


rm(data.path, treatments.path)


# convert specified columns to factors
data <- data %>%
  mutate(site = as.factor(site),
         layer = as.factor(layer),
         plot = as.factor(plot))

# restrict to times 0, 1, 4, since remaining times not yet measured;
# remove non-canopy photos
data <- data %>% filter(!TimeStep %in% c("2", "3"),
                        layer != "all")

summary(data) # note: 109 time 0, 109 time 1a+1b, 108 time 4; 109 L1, 109 L3, 108 L2
sort(data$plot[c(which(data$TimeStep==0))])
data[which(data$plot=="L3-5"),]
# photo tcDSCN0870 is an empty duplicate of L3-5
dim(data)
data %>% filter(photo != "tcDSCN0870") %>% dim() # cuts too many rows! should be 3*108 = 324 rows remaining
data %>% as.tbl() %>% filter(is.na(photo)) # the NAs do it
data %>% filter(!photo %in% c("tcDSCN0870")) %>% dim() # this works
data <- data %>% filter(!photo %in% c("tcDSCN0870"))

# still one extra row (should be 36 plot/site * 3 sites * 3 times = 324 total plots)
data %>%
  group_by(TimeStep, plot) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
# plot L1-2 is duplicated, for TimeStep 1(a)
data %>% filter(plot=="L1-2",
                TimeStep == "1a")
# prefer second pass, so remove the first hit to plot == "L1-2" & TimeStep == "1a"
dup <- which(data$plot == "L1-2" & data$TimeStep == "1a")[1]
data[dup,] # that's the one
# remove it
data <- data[-dup, ]
rm(dup)

head(data)
dim(data)

### prepare for analyses ####

#   consolidate certain taxa based on my ability to consistently ID them
#   simplify names 

# prepare to consolidate all brown encrusting algae: convert NA values to 0
data$`brown encrusting (red)`[which(is.na(data$`brown encrusting (red)`))]<-0
data$`Hildenbrandia`[which(is.na(data$`Hildenbrandia`))]<-0
data$`Petrocelis`[which(is.na(data$`Petrocelis`))]<-0
data$`Ralfsia`[which(is.na(data$`Ralfsia`))]<-0

# edit names
data <- data %>%
  rename(crust.red.brown = 'brown encrusting (red)',
         dead.barnacle = 'dead barnacle',
         littorine.sm = 'littorine - small (<1cm = <130px)',
         littorine.lg = 'littorine - large (>1cm = >130px)',
         crab.sm = 'crab - small',
         crab.lg = 'crab - large',
         epiphyte.red.brown.scuzzy = 'brown/red scuzzy epiphyte',
         hydroid = 'hydroid (Dynamena?)',
         loose.surface = 'sand/gravel/chips/debris')
colnames(data)

# consolidate encrusting algae, crabs, & ephemeral algae
data <- data %>%
  mutate(crust.red.brown = crust.red.brown + Hildenbrandia + Petrocelis + Ralfsia,
         crab = crab.lg + crab.sm,
         ephemeral = Ulva + Dumontia + epiphyte.red.brown.scuzzy + Porphyra) %>%
  select(-Hildenbrandia, -Petrocelis, -Ralfsia, -crab.lg, -crab.sm,
         -Ulva, -Dumontia, -epiphyte.red.brown.scuzzy, -Porphyra)

data$date <- as.Date(data$date, "%m/%d/%Y")

summary(data)

# clean up TimeStep, convert to factor
data$TimeStep[c(which(data$TimeStep=="1a"), which(data$TimeStep=="1b"))]<-"1"
data$TimeStep <- as.factor(data$TimeStep)
levels(data$TimeStep)

# more tidying up
data <- data %>%
  # remove uninformative columns (layer, photo, empty mastocarpus)
  # remove mobile inverts
  select(-layer, -photo, -Mastocarpus, 
         -littorine.lg, -littorine.sm, -crab, -Nucella) %>%
  # convert hydroid to dbl
  mutate(hydroid = as.numeric(hydroid))

summary(data)

### append plot treatments to data ####

head(treatments)

# change "removal" column to "diversity(sp_removed)' for now, for clarity
# also change all column names to lower case
treatments <- treatments %>% 
  rename(removal = 'diversity(sp_removed)') %>%
  rename_all(tolower)

# rename the cage-ctrl treatment as "roof"
treatments$pred[which(treatments$pred == "cage-ctrl")] <- "roof"

# convert to factor
treatments <- treatments %>%
  mutate(site = as.factor(site),
         id = as.factor(id),
         removal = as.factor(removal),
         pred = as.factor(pred))

summary(treatments)

# looks as it should.
# rename id column as "plot"
treatments <- treatments %>% rename(plot = id)

# join the main data frame
data <- data %>% left_join(treatments, by = c("site" = "site", "plot" = "plot"))
rm(treatments)

summary(data)
str(data)

# arrange by plot and time
data <- data %>%
  arrange(plot, TimeStep)


### create data frame showing cover of each target species over time
# melt
mdat <- melt(data, id.vars= c("site","plot","TimeStep","date","pred","removal"),
           measure.vars = c("Fucus","Mytilus","barnacle"), 
           variable.name="species", value.name="pct.cover")
head(mdat)
summary(mdat)
str(mdat)

mdat$sp.T<-paste(mdat$species,mdat$TimeStep,sep="-")
mdat$sp.T<-factor(mdat$sp.T)

# recast
dat<-dcast(mdat,site+plot+pred+removal~sp.T,value.var="pct.cover")
summary(dat)

which(dat$removal=="A"|dat$removal=="B")

#### test reinvasion: is end cover greater than immediate post-removal cover? ####

## Fucus
# proportion of plots reinvaded where Fucus was target removal
sum((dat$`Fucus-4`[which(dat$removal=="A")]-dat$`Fucus-1`[which(dat$removal=="A")])>0)/27
# proportion of plots reinvaded where barnacles or Mytilus were targets of removal
sum((dat$`Fucus-4`[which(dat$removal=="B"|dat$removal=="M")]-
       dat$`Fucus-1`[which(dat$removal=="B"|dat$removal=="M")])>0)/54
# proportion of plots reinvaded where nothing was removed 
sum((dat$`Fucus-4`[which(dat$removal=="none")]-dat$`Fucus-1`[which(dat$removal=="none")])>0)/27

## barnacle
# proportion of plots reinvaded where barnacles were target removal
sum((dat$`barnacle-4`[which(dat$removal=="B")]-dat$`barnacle-1`[which(dat$removal=="B")])>0)/27
# proportion of plots reinvaded where Fucus or Mytilus were targets of removal
sum((dat$`barnacle-4`[which(dat$removal=="A"|dat$removal=="M")]-
       dat$`barnacle-1`[which(dat$removal=="A"|dat$removal=="M")])>0)/54
# proportion of plots reinvaded where nothing was removed 
sum((dat$`barnacle-4`[which(dat$removal=="none")]-dat$`barnacle-1`[which(dat$removal=="none")])>0)/27

## Mytilus
# proportion of plots reinvaded where Mytilus was target removal
sum((dat$`Mytilus-4`[which(dat$removal=="M")]-dat$`Mytilus-1`[which(dat$removal=="M")])>0)/27
# proportion of plots reinvaded where Fucus or barnacles were targets of removal
sum((dat$`Mytilus-4`[which(dat$removal=="B"|dat$removal=="A")]-
       dat$`Mytilus-1`[which(dat$removal=="B"|dat$removal=="A")])>0)/54
# proportion of plots reinvaded where nothing was removed 
sum((dat$`Mytilus-4`[which(dat$removal=="none")]-dat$`Mytilus-1`[which(dat$removal=="none")])>0)/27


# add regrowth amount of each species (from time 1 to time 4) to the data frame
dat$re.A<-1*(dat$`Fucus-4`-dat$`Fucus-1`>0)
dat$re.B<-1*(dat$`barnacle-4`-dat$`barnacle-1`>0)
dat$re.M<-1*(dat$`Mytilus-4`-dat$`Mytilus-1`>0)

# what was unintentional Fucus removal from barnacle plots?
min(dat$`Fucus-1`[which(dat$removal=="B")]-dat$`Fucus-0`[which(dat$removal=="B")])
max(dat$`Fucus-1`[which(dat$removal=="B")]-dat$`Fucus-0`[which(dat$removal=="B")])

mean((dat$`Fucus-1`[which(dat$removal=="B")]-dat$`Fucus-0`[which(dat$removal=="B")])/
       dat$`Fucus-0`[which(dat$removal=="B")])*100
min((dat$`Fucus-1`[which(dat$removal=="B")]-dat$`Fucus-0`[which(dat$removal=="B")])/
       dat$`Fucus-0`[which(dat$removal=="B")])*100
max((dat$`Fucus-1`[which(dat$removal=="B")]-dat$`Fucus-0`[which(dat$removal=="B")])/
       dat$`Fucus-0`[which(dat$removal=="B")])*100

t.test(dat$`Fucus-0`[which(dat$removal=="B")],dat$`Fucus-1`[which(dat$removal=="B")],
       alternative="two.sided",paired=TRUE)


t.test(dat$`Fucus-4`[which(dat$removal=="none")],dat$`Fucus-1`[which(dat$removal=="none")],
       alternative="two.sided",paired=TRUE)

t.test(dat$`Fucus-4`[which(dat$removal=="none" & dat$pred=="open")],
       dat$`Fucus-1`[which(dat$removal=="none" & dat$pred=="open")],
       alternative="two.sided",paired=TRUE)


a.growth <- dat$`Fucus-4`[which(dat$removal=="none")] - dat$`Fucus-1`[which(dat$removal=="none")]
preds <- dat$pred[which(dat$removal=="none")]
summary(aov(a.growth ~ preds))

table(dat$removal, dat$re.A)
re.A.vec = tapply(dat$re.A,dat$removal,sum)
tot.A.vec = tapply(1+dat$re.A-dat$re.A,dat$removal,sum)
prop.test(re.A.vec, tot.A.vec)

table(dat$removal, dat$re.B)
re.B.vec = tapply(dat$re.B,dat$removal,sum)
tot.B.vec = tapply(1+dat$re.B-dat$re.B,dat$removal,sum)
prop.test(re.B.vec, tot.B.vec)

table(dat$removal, dat$re.M)
re.M.vec = tapply(dat$re.M,dat$removal,sum)
tot.M.vec = tapply(1+dat$re.M-dat$re.M,dat$removal,sum)
prop.test(re.M.vec, tot.M.vec)


binom.confint(x = 23, n = 27, methods="wilson", conf.level = 0.95) #A
binom.confint(x = 27, n = 27, methods="wilson", conf.level = 0.95) #B
binom.confint(x = 26, n = 27, methods="wilson", conf.level = 0.95) #M

table(dat$pred, dat$re.A) 
table(dat$pred, dat$re.B) 
table(dat$pred, dat$re.M)

prop.test(23,27, p=0.5, alternative="greater",correct = TRUE)
prop.test(27,27, p=0.5, alternative="greater", correct=TRUE)
prop.test(26,27, p=0.5, alternative="greater", correct=TRUE)

props <- numeric(length = 3)
lower <- numeric(length = 3)
upper <- numeric(length = 3)
df <- data.frame(cbind(props, lower, upper))
df$successes <- c(23,27,26)
df$props <- df$successes/27
for(i in 1:3){
  df$lower[i] <- prop.test(df$successes[i],27, p=0.5, alternative="greater", correct=TRUE)$conf.int[1]
}  
for(i in 1:3){
  df$upper[i] <- prop.test(df$successes[i],27, p=0.5, alternative="greater", correct=TRUE)$conf.int[2]
}  
df$species<-factor(c("Fucus", "Barnacle","Mussel"), levels=c("Fucus", "Barnacle", "Mussel"))


ggplot(data=df,mapping=aes(y=props,x=species))+
  geom_point()+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.1)+
  # ggtitle("Regrowth of focal species following removal")+
  labs(y="Proportion of plots with regrowth",x="Focal species")+
  ylim(0,1)+
  geom_hline(yintercept = 0, linetype=2)+
  theme_bw()
  # theme(axis.line = element_line(colour = "black"),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       panel.border = element_blank(),
  #       panel.background = element_blank(),
  #       axis.ticks=element_blank(),
  #       axis.line.x =element_blank()) 

ps <- c(
  prop.test(23,27, p=0.5, alternative="greater", correct=TRUE)$p.value,
  prop.test(27,27, p=0.5, alternative="greater", correct=TRUE)$p.value,
  prop.test(26,27, p=0.5, alternative="greater", correct=TRUE)$p.value
)

p.adjust(ps, method="holm")

dat.A <- dat[which(dat$removal=="A"),]
dat.B <- dat[which(dat$removal=="B"),]
dat.M <- dat[which(dat$removal=="M"),]

table(dat.A$pred, dat.A$re.A) # reinvasions in 8-8-7 for cage-open-roof
table(dat.B$pred, dat.B$re.B) # reinvasions in 9-9-9 for cage-open-roof
table(dat.M$pred, dat.M$re.M) # reinvasions in 9-9-8 for cage-open-roof

prop.test(c(8,8,7),c(9,9,9))
prop.test(c(9,9,8),c(9,9,9))
prop.test(c(9,9,9),c(9,9,9))

# regrowth amount - arithmetic (linear) growth
dat$d.A <- dat$`Fucus-4`-dat$`Fucus-1`
dat$d.B <- dat$`barnacle-4`-dat$`barnacle-1`
dat$d.M <- dat$`Mytilus-4`-dat$`Mytilus-1`
rA <- which(dat$removal=="A")
rB <- which(dat$removal=="B")
rM <- which(dat$removal=="M")
rN <- which(dat$removal == "none")

mean(dat$d.A[rA])
mean(dat$d.B[rB])
mean(dat$d.M[rM])
sd(dat$d.A[rA])/sqrt(27)
sd(dat$d.B[rB])/sqrt(27)
sd(dat$d.M[rM])/sqrt(27)

regrowth <- numeric(length=3)
df2 <- data.frame(regrowth)
df2$regrowth <- c(mean(dat$d.A[rA]),
                mean(dat$d.B[rB]),
                mean(dat$d.M[rM]))
df2$se <- c(sd(dat$d.A[rA]) / sqrt(27),
          sd(dat$d.B[rB]) / sqrt(27),
          sd(dat$d.M[rM]) / sqrt(27))
df2$species <- df$species

ggplot(df2, aes(y=regrowth, x=species))+
  geom_point()+
  geom_errorbar(aes(ymin=regrowth-se, ymax=regrowth+se), width=.1)+
  # ggtitle("Regrowth of focal species following removal")+
  labs(y="Increase in cover after removal \n (percentage points)",x="Focal species")+
  geom_hline(yintercept = 0, linetype=2)+
  theme_bw()

df3 <- data.frame(matrix(nrow=6))
df3$growth <- c(mean(dat$d.A[rA]),
              mean(dat$d.B[rB]),
              mean(dat$d.M[rM]),
              mean(dat$d.A[-rA]),
              mean(dat$d.B[-rB]),
              mean(dat$d.M[-rM]))
df3$se <- c(sd(dat$d.A[rA])/ sqrt(81),
          sd(dat$d.B[rB]) / sqrt(81),
          sd(dat$d.M[rM]) / sqrt(81),
          sd(dat$d.A[-rA]) / sqrt(81),
          sd(dat$d.B[-rB]) / sqrt(81),
          sd(dat$d.M[-rM]) / sqrt(81))

df3[,1] <- NULL
df3$species <- rep(df2$species, 2)
df3$Status <- factor(c(rep("invader", 3), rep("resident", 3)))
df3

ggplot(df3, aes(y=growth, x=species, col=Status))+
  geom_point()+
  geom_errorbar(aes(ymin=growth-se, ymax=growth+se), width=.1)+
  # ggtitle("Regrowth of focal species following removal")+
  labs(y="Change in cover after treatment \n (percentage points)",x="Focal species")+
  geom_hline(yintercept = 0, linetype=2)+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks=element_blank(),
        axis.line.x =element_blank(),
        text = element_text(size=16)) 

df4 <- data.frame(matrix(nrow=12))
df4$growth <- c(mean(dat$d.A[rA]),
              mean(dat$d.B[rA]),
              mean(dat$d.M[rA]),
              mean(dat$d.A[rB]),
              mean(dat$d.B[rB]),
              mean(dat$d.M[rB]),
              mean(dat$d.A[rM]),
              mean(dat$d.B[rM]),
              mean(dat$d.M[rM]),
              mean(dat$d.A[rN]),
              mean(dat$d.B[rN]),
              mean(dat$d.M[rN]))
df4$se <- c(sd(dat$d.A[rA])/sqrt(27),
          sd(dat$d.B[rA])/sqrt(27),
          sd(dat$d.M[rA])/sqrt(27),
          sd(dat$d.A[rB])/sqrt(27),
          sd(dat$d.B[rB])/sqrt(27),
          sd(dat$d.M[rB])/sqrt(27),
          sd(dat$d.A[rM])/sqrt(27),
          sd(dat$d.B[rM])/sqrt(27),
          sd(dat$d.M[rM])/sqrt(27),
          sd(dat$d.A[rN])/sqrt(27),
          sd(dat$d.B[rN])/sqrt(27),
          sd(dat$d.M[rN])/sqrt(27))

df4[,1] <- NULL
df4$species <- rep(df2$species, 4)
df4$Invader <- factor(c(rep("Fucus", 3), rep("Barnacle", 3), rep("Mussel", 3), rep("Control",3)))
df4

dodge = position_dodge(width=0.2)
ggplot(df4, aes(y=growth, x=species, col=Invader))+
  geom_point(position=dodge)+
  geom_errorbar(aes(ymin=growth-se, ymax=growth+se), width=.1, position=dodge)+
  # ggtitle("Regrowth of focal species following removal")+
  labs(y="Change in cover after treatment \n (percentage points)",x="Focal species")+
  geom_hline(yintercept = 0, linetype=2)+
  scale_color_manual(values = c("#F8766D", "#C77CFF", "#7CAE00", "#00BFC4"))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks=element_blank(),
        axis.line.x =element_blank(),
        text = element_text(size=16)) 


# dat$r.A<-
  ((dat$`Fucus-4`+5)/(dat$`Fucus-1`+5))


### Ordination ####

### nMDS

names(data)
community <- data[,6:16]
id.tags <- data[,c(1:4,17:18)]


community.mds<-metaMDS(community,distance="bray",trymax=50,k=3)
### check stress plot
stressplot(community.mds)

# community.mds.j<-metaMDS(community,distance="jaccard",binary=TRUE, trymax=50,k=3)
# stressplot(community.mds.j)
# 
# pro<-procrustes(community.mds,community.mds.j)
# pro
# plot(pro)
# plot(pro,kind=2)

### visualize plot
plot(community.mds)
treat=id.tags$removal
ordiplot(community.mds,type="n")
ordihull(community.mds,groups=treat,draw="polygon",col="grey90",label=F)
# orditorp(example_NMDS,display="species",col="red",air=0.01)
# orditorp(example_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
#          air=0.01,cex=1.25)


t0<-seq(from=1,by=3,length=108)
t1<-seq(from=2,by=3,length=108)
t4<-seq(from=3,by=3,length=108)



# community0.mds2<-metaMDS(data[t0,7:20],k=2)
# stressplot(community0.mds2)
community0.mds3<-metaMDS(community[t0,],k=3,trymax=100)
# comm.dist0<-vegdist(data[t0,6:20])
stressplot(community0.mds3)
# stressplot(community0.mds3,comm.dist0)

# community1.mds2<-metaMDS(data[t1,7:20],k=2)
# stressplot(community1.mds2)
community1.mds3<-metaMDS(community[t1,],k=3)
stressplot(community1.mds3)

# community4.mds2<-metaMDS(data[t4,7:20],k=2,trymax=100)
# stressplot(community4.mds2)
community4.mds3<-metaMDS(community[t4,],k=3, trymax = 100)
stressplot(community4.mds3)


treat<-id.tags$removal[t0]
colors<-c(rep(c("gray90","light blue","tomato"),10),
          rep(c("gray60","cornflower blue","orangered"),10),
          rep(c("gray40","blue","red3"),10),
          rep(c("black","dark blue","red4"),10))
colvec <- c("green4", "red2", "mediumblue","orange")
pchs<-c(2,8,16)
gr.removal<-factor(data$removal)
gr.time<-factor(data$TimeStep)

plot(community.mds, type = "n", display = "sites")
orditorp(community.mds,display="species",col="black",air=0.01)
points(community.mds, display = "sites", pch = pchs[gr.time], col = colvec[gr.removal])
ordihull(community.mds, display = "sites", groups=gr.removal, col = colvec)
ordiellipse(community.mds, display = "sites", groups=gr.removal, col = colvec, kind="sd")
legend("topright", legend=levels(gr.time), bty = "n", col= c("black"), pch = pchs,title="TimeStep")
legend("bottomright", legend = levels(gr.removal), bty = "n", col = colvec, pch=c(20),title="Removal")

ordiplot(community0.mds3,type="n",main="Initial")
# ordihull(community0.mds3,groups=treat,col= colvec,label=F)
ordiellipse(community0.mds3, display = "sites", groups=treat, col = colvec, kind="sd",lty=2)
ordiellipse(community0.mds3, display = "sites", groups=treat, col = colvec, kind="se")
orditorp(community0.mds3,display="species",col="red",air=0.01)
# orditorp(community0.mds3,display="sites",col=colvec[treat],air=0.01)
legend("topright", legend = levels(gr.removal), bty = "n", col = colvec, pch=c(20),title="Removal")

ordiplot(community1.mds3,type="n",main="Treatment")
# ordihull(community1.mds3,groups=treat,col= colvec,label=F)
ordiellipse(community1.mds3, display = "sites", groups=treat, col = colvec, kind="sd",lty=2)
ordiellipse(community1.mds3, display = "sites", groups=treat, col = colvec, kind="se")
orditorp(community1.mds3,display="species",col="red",air=0.01)
# orditorp(community1.mds3,display="sites",col=colvec[treat],air=0.01)
ordispider(community1.mds3, display = "sites", groups=treat, col = colvec, spiders="centroid")
legend("bottomright", legend = levels(gr.removal), bty = "n", col = colvec, pch=c(20),title="Removal")

ordiplot(community4.mds3,type="n",main="Recovery")
# ordihull(community4.mds3,groups=treat,col= colvec,label=F)
ordiellipse(community4.mds3, display = "sites", groups=treat, col = colvec, kind="sd",lty=2)
ordiellipse(community4.mds3, display = "sites", groups=treat, col = colvec, kind="se")
ordispider(community4.mds3, display = "sites", groups=treat, col = colvec, spiders="centroid")
orditorp(community4.mds3,display="species",col="red",air=0.01)
# orditorp(community4.mds3,display="sites",col=colvec[treat],air=0.01)
ordibar(community4.mds3, groups=treat, display = "species", kind = "sd", col = colvec, label = FALSE)
legend("bottomright", legend = levels(gr.removal), bty = "n", col = colvec, pch=c(20),title="Removal")

# library(vegan3d)
# ordiplot3d(community0.mds3,choices=c(1:3),col=colvec)
# ordiplot3d(community1.mds3,choices=c(1:3),col=colvec)
# ordiplot3d(community4.mds3,choices=c(1:3),col=colvec)
# orglpoints(community0.mds3,display="sites",col=colvec,air=0.01,cex=1.25)
# orglpoints(community1.mds3,display="sites",col=colvec,air=0.01,cex=1.25)
# orglpoints(community4.mds3,display="sites",col=colvec,air=0.01,cex=1.25)

MDS1<-community.mds$points[,1]
MDS2<-community.mds$points[,2]
MDS3<-community.mds$points[,3]

NMDS <- data.frame(MDS1 = MDS1, MDS2 = MDS2, MDS3 = MDS3, 
                   Plot = data$plot, Site = data$site, TimeStep = data$TimeStep,
                   Removal = data$removal, Pred = data$pred)
summary(NMDS)
NMDS$Removal<-as.character(NMDS$Removal)
NMDS$Removal[which(NMDS$Removal=="A")]<-"Fucus"
NMDS$Removal[which(NMDS$Removal=="B")]<-"Barnacle"
NMDS$Removal[which(NMDS$Removal=="M")]<-"Mussel"
NMDS$Removal[which(NMDS$Removal=="none")]<-"No Removal"
NMDS$Removal<-factor(NMDS$Removal, levels=c("Fucus","Barnacle","Mussel","No Removal"))

centroids <- aggregate(cbind(MDS1,MDS2,MDS3)~Removal*Pred*TimeStep,NMDS,mean)
centroids.r <- aggregate(cbind(MDS1,MDS2,MDS3)~Removal*TimeStep,NMDS,mean)

NMDS.species<-data.frame(scores(community.mds, display="species"))
rownames(NMDS.species)[5]<-"encrusting"
rownames(NMDS.species)[7]<-"debris"


# Site groupings
ggplot(NMDS[t0,], aes(x=MDS1, y=MDS2, col=Site)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "Initial (pre-treatment)")+
  coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))+
  geom_text(inherit.aes=FALSE, data=NMDS.species, label=rownames(NMDS.species), 
            aes(x=NMDS1, y=NMDS2), color="red")

ggplot(NMDS[t1,], aes(x=MDS1, y=MDS2, col=Site)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "Treatment")+
  coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))+
  geom_text(inherit.aes=FALSE, data=NMDS.species, label=rownames(NMDS.species), 
            aes(x=NMDS1, y=NMDS2), color="red")

ggplot(NMDS[t4,], aes(x=MDS1, y=MDS2, col=Site)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "Recovery")+
  coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))+
  geom_text(inherit.aes=FALSE, data=NMDS.species, label=rownames(NMDS.species), 
            aes(x=NMDS1, y=NMDS2), color="red")

# Diversity groupings
ggplot(NMDS[t0,], aes(x=MDS1, y=MDS2, col=Removal)) +
  geom_point() +
  stat_ellipse() +
  scale_color_manual(values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF"))+
  theme_bw() +
  labs(title = "Initial (pre-treatment)")+
  coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))+
  geom_text(inherit.aes=FALSE, data=NMDS.species[c(1,2,4),], label=rownames(NMDS.species)[c(1,2,4)], 
            aes(x=NMDS1, y=NMDS2), color="red")

ggplot(NMDS[t1,], aes(x=MDS1, y=MDS2, col=Removal)) +
  geom_point() +
  scale_color_manual(values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF"))+
  stat_ellipse(level=.95) +
  theme_bw() +
  labs(title = "Treatment")+
  coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))+
  geom_text(inherit.aes=FALSE, data=NMDS.species[c(1,2,4),], label=rownames(NMDS.species)[c(1,2,4)], 
            aes(x=NMDS1, y=NMDS2), color="red")

ggplot(NMDS[t4,], aes(x=MDS1, y=MDS2, col=Removal)) +
  geom_point() +
  scale_color_manual(values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF"))+
  stat_ellipse() +
  theme_bw() +
  # geom_errorbar(aes(ymin=lower, ymax=upper))+
  labs(title = "Recovery")+
  coord_cartesian(xlim = c(-1,1), ylim = c(-1,1))+
  geom_text(inherit.aes=FALSE, data=NMDS.species[c(1,2,4),], label=rownames(NMDS.species)[c(1,2,4)], 
            aes(x=NMDS1, y=NMDS2), color="red")


# Predator groupings
ggplot(NMDS[t0,], aes(x=MDS1, y=MDS2, col=Pred)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "Initial (pre-treatment)")
ggplot(NMDS[t1,], aes(x=MDS1, y=MDS2, col=Pred)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "Initial (pre-treatment)")
ggplot(NMDS[t4,], aes(x=MDS1, y=MDS2, col=Pred)) +
  geom_point() +
  stat_ellipse() +
  theme_bw() +
  labs(title = "Initial (pre-treatment)")

# # Pred*Diversity groupings
# ggplot(NMDS[t0,], aes(x=MDS1, y=MDS2, col=interaction(Removal,Pred))) +
#   geom_point() +
#   stat_ellipse() +
#   theme_bw() +
#   labs(title = "Initial (pre-treatment)")
# ggplot(NMDS[t0,], aes(x=MDS3, y=MDS2, col=interaction(Removal,Pred))) +
#   geom_point() +
#   stat_ellipse() +
#   theme_bw() +
#   labs(title = "Initial (pre-treatment)")
# ggplot(NMDS[t0,], aes(x=MDS1, y=MDS3, col=interaction(Removal,Pred))) +
#   geom_point() +
#   stat_ellipse() +
#   theme_bw() +
#   labs(title = "Initial (pre-treatment)")

ggplot(NMDS[t0,], aes(x=MDS1, y=MDS2, col=Removal,pch=Pred)) +
  # geom_point() +
  scale_color_manual(values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF"))+
  scale_linetype_manual(values = c(1,3,2))+
  stat_ellipse(aes(linetype=Pred)) +
  theme_bw() +
  labs(title = "Initial (pre-treatment)")+
  geom_point(data=centroids[which(centroids$TimeStep==0),],size=2)+
  coord_cartesian(xlim = c(-1,1.1), ylim = c(-1.1,1.5))+
  geom_text(inherit.aes=FALSE, data=NMDS.species, label=rownames(NMDS.species), 
            aes(x=NMDS1, y=NMDS2), color="red")

ggplot(NMDS[t1,], aes(x=MDS1, y=MDS2, col=Removal,pch=Pred)) +
  # geom_point() +
  scale_color_manual(values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF"))+
  scale_linetype_manual(values = c(1,3,2))+
  stat_ellipse(aes(linetype=Pred)) +
  theme_bw() +
  labs(title = "Treatment")+
  geom_point(data=centroids[which(centroids$TimeStep==1),],size=2)+
  coord_cartesian(xlim = c(-1,1.1), ylim = c(-1.1,1.5))+
  geom_text(inherit.aes=FALSE, data=NMDS.species, label=rownames(NMDS.species), 
            aes(x=NMDS1, y=NMDS2), color="red")

ggplot(NMDS[t4,], aes(x=MDS1, y=MDS2, col=Removal,pch=Pred)) +
  # geom_point() +
  scale_color_manual(values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF"))+
  scale_linetype_manual(values = c(1,3,2))+
  stat_ellipse(aes(linetype=Pred)) +
  theme_bw() +
  labs(title = "Recovery")+
  geom_point(data=centroids[which(centroids$TimeStep==4),],size=2, position="jitter")+
  coord_cartesian(xlim = c(-1,1.1), ylim = c(-1.1,1.5))+
  geom_text(inherit.aes=FALSE, data=NMDS.species, label=rownames(NMDS.species), 
            aes(x=NMDS1, y=NMDS2), color="red")


# figure out ggplot default colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)
plot(1:n, pch = 16, cex = 2, col = cols)

# figure out ggplot default line types
lines <- data.frame(x = 1:2, y = rep(13:1, each = 2), grp = factor(rep(1:13, each = 2)))
lines
# plot
ggplot(data = lines, aes(x = x, y = y, linetype = grp)) +
  geom_line() +
  geom_text(aes(x = 0.95, label = grp)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none")

community.mds
a0<-anosim(data[t0,6:16],grouping=data$removal[t0])
a1<-anosim(data[t1,6:16],grouping=data$removal[t1])
a4<-anosim(data[t4,6:16],grouping=data$removal[t4])
summary(a0)
summary(a1)
summary(a4)
plot(a0)
plot(a1)
plot(a4)

### PERMANOVA

# full model
adonis.full<-adonis(vegdist(data[,c(5:16)],method="bray")
       ~TimeStep*pred*removal, strata=data$site, data=data, permutations=999)

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

interaction.plot(data.0$pred, data.0$removal, data.0$Fucus,ylab="Fucus cover",xlab="predator treatment")
interaction.plot(data.0$pred, data.0$removal, data.0$barnacle,ylab="barnacle cover",xlab="predator treatment")
interaction.plot(data.0$pred, data.0$removal, data.0$Mytilus,ylab="Mytilus cover",xlab="predator treatment")
interaction.plot(data.0$pred, data.0$removal, data.0$crust.red.brown,ylab="encrusting",xlab="predator treatment")

tapply(data.0$Fucus,data.0$pred:data.0$removal,mean) # noisy. greatest cover in roof:none, lowest in roof:mussel
tapply(data.0$barnacle,data.0$pred:data.0$removal,mean) # barnacle cover lowest in cage:none, comparable elsewhere
tapply(data.0$Mytilus,data.0$pred:data.0$removal,mean) # mussel cover high in (CAGE & open):NONE vs (CAGE & open):rest and low in ROOF:NONE vs roof:rest

mean(data[which(TimeStep == 0)]$bare)
mean(data[which(TimeStep == 0 & removal == "none")]$bare)
mean(data[which(TimeStep == 0 & removal != "none")]$bare)
mean(data[which(TimeStep == 1 & removal != "none")]$bare)
mean(data[which(TimeStep == 4 & removal == "none")]$bare)


# mixed effects modeling?

m1 <- lmer(Mytilus~removal * pred + (1|site) + (1|TimeStep), data=data)
m2 <- lmer(Mytilus~removal + pred + (1|site) + (1|TimeStep), data=data)
KRmodcomp(m1, m2)

b1 <- lmer(barnacle~removal * pred + (1|site) + (1|TimeStep), data=data)
b2 <- lmer(barnacle~removal + pred + (1|site) + (1|TimeStep), data=data)
KRmodcomp(b1, b2)

f1 <- lmer(Fucus~removal * pred + (1|site) + (1|TimeStep), data=data)
f2 <- lmer(Fucus~removal + pred + (1|site) + (1|TimeStep), data=data)
KRmodcomp(f1, f2)

species_names = c("barnacle" = "Barnacle", 
                  "Fucus" = "Fucus", 
                  "Mytilus" = "Mytilus", 
                  "Control" = "Control"
                  )
dodge=position_dodge(width=.1)

ggplot(mdat, aes(x=TimeStep,y=pct.cover,color=removal, pch=pred))+
  stat_summary(fun.y = mean,geom="point",size=2,position=dodge)+
  stat_summary(fun.y = mean,geom="line",aes(group=interaction(removal, pred)),linetype=2,position=dodge)+
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.1,position=dodge)+
  facet_wrap(~species, nrow=3, labeller = as_labeller(species_names))+
  scale_x_discrete(name = "Time Step", labels = c("Initial", "Treatment", "Recovery"))+
  ylab("Percent Cover")+
  scale_color_manual(name="Removal",labels=c("Fucus","Barnacle", "Mytilus", "Control"),
                     values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF"))+
  theme_bw(base_size=18)

mdat2<-melt(data, id.vars= c("site","plot","TimeStep","date","pred","removal"),
           measure.vars = "bare", 
           variable.name="species", value.name="pct.cover")
ggplot(mdat2, aes(x=TimeStep, y=pct.cover, color=removal))+
  stat_summary(fun.y=mean, geom="point", size=2, position=dodge)+
  stat_summary(fun.y = mean,geom="line",aes(group=removal),linetype=2,position=dodge)+
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.1,position=dodge)+
  xlab("Time Step")+
  ylab("Percent Bare Space")+
  scale_color_manual(name="Removal",labels=c("Fucus","Barnacle", "Mytilus", "Control"),
                     values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF"))+
  scale_x_discrete(labels=c("Initial", "Treatment", "Recovery"))+
  theme_bw(base_size=18)+
  theme(axis.text.x = element_text(hjust=0.5, vjust=1))

dodge=position_dodge(width=.1)
ggplot(mdat2, aes(x=TimeStep, y=pct.cover, color=removal, pch=pred))+
  stat_summary(fun.y=mean, geom="point", size=2, position=dodge, aes(group=interaction(removal, pred)))+
  stat_summary(fun.y = mean,geom="line",aes(group=interaction(removal, pred)),linetype=2,position=dodge)+
  stat_summary(fun.data = mean_se,geom="errorbar",width=0.1,position=dodge, aes(group=interaction(removal, pred)))+
  xlab("Time Step")+
  ylab("Percent Bare Space")+
  scale_color_manual(name="Removal",labels=c("Fucus","Barnacle", "Mytilus", "Control"),
                     values = c("#7CAE00", "#F8766D", "#00BFC4", "#C77CFF"))+
  scale_x_discrete(labels=c("Initial", "Treatment", "Recovery"))+
  theme_bw(base_size=18)+
  theme(axis.text.x = element_text(hjust=0.5, vjust=1))



data %>%
  filter(TimeStep == 0) %>%
  group_by(removal, pred) %>%
  # filter(removal == "none" & pred == "roof") %>%
  select(Fucus, barnacle, Mytilus) %>%
  summarise(a_avg = mean(Fucus),
            b_avg = mean(barnacle),
            m_avg = mean(Mytilus),
            inverts = b_avg + m_avg) %>%
  arrange(inverts)

ds<-data.0[-which(data.0$removal=="none" & data.0$pred=="roof"),]
adonis(vegdist(ds[,c(6:16)],method="bray")
       ~pred*removal, strata=ds$site, data=ds, permutations=999) # sig/ marginally sig


data %>%
  filter(TimeStep == 0) %>%
  group_by(removal, pred) %>%
  # filter(removal == "none" & pred == "roof") %>%
  select(Fucus, barnacle, Mytilus) %>%
  filter(removal == "A" & pred == "open")

data %>%
  filter(TimeStep == 0) %>%
  group_by(removal, pred) %>%
  # filter(removal == "none" & pred == "roof") %>%
  select(Fucus, barnacle, Mytilus) %>%
  filter(!(removal == "A" & pred == "open"))

data[which(TimeStep == 0 & removal == "A" & pred == "open")]
data[which(TimeStep == 0 & removal == "none" & pred == "roof")]
