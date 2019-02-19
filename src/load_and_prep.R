library(tidyverse)
library(lubridate)

### load data #####################

# data import wrapper fxn with "...", just for fun!
path.read.csv <- function(...){
  data_path <- file.path(...)
  df <- read_csv(data_path)
  return(df)
}

data <- path.read.csv("data", "counts.csv")

dim(data)
str(data)
head(data)
tail(data)
# there are some extra rows, will deal with them soon

## load treatments for each plot
treatments <- path.read.csv("data", "treatments.csv")
dim(treatments)

head(treatments)
tail(treatments, n = 15)
# retain Site & ID as key columns, diversity and pred treatment columns; only first 108 rows have useful data
treatments <- treatments[1:108,c(1,3:5)]
treatments
tail(treatments)

###### address the extra rows ##################

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
data %>% filter(is.na(photo)) # the NAs do it -- problems w/ comparing NA to a value!!
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

####### prepare for analyses ###################

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

#   consolidate certain taxa based on my ability to consistently ID them

# prepare to consolidate all brown encrusting algae: convert NA values to 0
data$crust.red.brown[which(is.na(data$crust.red.brown))] <- 0
data$Hildenbrandia[which(is.na(data$Hildenbrandia))] <- 0
data$Petrocelis[which(is.na(data$Petrocelis))] <- 0
data$Ralfsia[which(is.na(data$Ralfsia))] <- 0

# consolidate encrusting algae, crabs, & ephemeral algae
data <- data %>%
  mutate(crust.red.brown = crust.red.brown + Hildenbrandia + Petrocelis + Ralfsia,
         crab = crab.lg + crab.sm,
         ephemeral = Ulva + Dumontia + epiphyte.red.brown.scuzzy + Porphyra) %>%
  select(-Hildenbrandia, -Petrocelis, -Ralfsia, -crab.lg, -crab.sm,
         -Ulva, -Dumontia, -epiphyte.red.brown.scuzzy, -Porphyra)

# consolidate TimeStep 1a and 1b as "1"
data$TimeStep[c(which(data$TimeStep=="1a"), which(data$TimeStep=="1b"))]<-"1"

## type conversions: factors, dates
data <- data %>%
  mutate(site = as.factor(site),
         layer = as.factor(layer),
         plot = as.factor(plot),
         TimeStep = as.factor(TimeStep),
         date = mdy(date))

summary(data)

# more tidying up
data <- data %>%
  # remove uninformative columns (layer, photo, empty mastocarpus)
  # remove mobile inverts
  select(-layer, -photo, -Mastocarpus, 
         -littorine.lg, -littorine.sm, -crab, -Nucella) %>%
  # convert hydroid to dbl
  mutate(hydroid = as.numeric(hydroid))

summary(data)

##### join plot treatments to data ###############

head(treatments)

# edit column names
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

str(data)

# arrange by plot and time
data <- data %>%
  arrange(plot, TimeStep)

summary(data)
str(data)
