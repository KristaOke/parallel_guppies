# GLMMs for parallel_guppies! 
# 2020-08-31 AH

# 2020-09-08 ah
# We need to collapse categories - TraitType1/TraitType2
# Probably should exclude common garden (check cgtally)
# Maybe include morphometric for CG?

# For collapsing traits... 
# Diet + Other + Physiology + Life History
# But for "both" physiology makes up a bigger proportion in the North
# For "both" maybe we can only look at South

# 2020-09-09 ah
# Collapsed the traits and re-ran the models... 

# 2020-09-14 ah
# Is sex.mod6 the best model? 
# Deleted all models with 'Paired' as an effect

# changed file names -- NEW FILES ON DRIVE

## 2021-04-03 AH - I've added a bunch and shoved the old stuff the the bottom

# LIBRARIES ---- 
library(plyr)
library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)
library(bbmle)
library(gridExtra)
library(car)
library(corrplot)
library(tidyverse)
library(fitdistrplus)
library(tidyverse)
library(sjPlot)

# set working directory
# setwd("")
# wd<- getwd()


## import and tidy 

# These are all the updated sheets on the Drive
# spreadsheet.data is the data extracted for the meta-analysis
# R2.data.among are the output of the ANOVA loops

spreadsheet.data <- read.csv(paste(wd,'/Data/MetaData.csv',sep=""), header=TRUE, sep=",")
R2.data.among <- read.csv(paste(wd,'/Data/TraitR2_among.csv',sep=""), header=TRUE, sep=",")
R2.data.south <- read.csv(paste(wd,'/Data/TraitR2_south.csv',sep=""), header=TRUE, sep=",")
R2.data.intro <- read.csv(paste(wd,'/Data/TraitR2_intro.csv',sep=""), header=TRUE, sep=",")
R2.data.caroni <- read.csv(paste(wd,'/Data/TraitR2_Caroni.csv',sep=""), header=TRUE, sep=",")
R2.data.among.drainage <- read.csv(paste(wd,'/Data/TraitR2_Among_Drainage.csv',sep=""), header=TRUE, sep=",")

str(spreadsheet.data)
str(R2.data.among)

spreadsheet.data$StudyID <- as.factor(spreadsheet.data$StudyID)
spreadsheet.data$Collection_start <- as.factor(spreadsheet.data$Collection_start)
spreadsheet.data$Collection_end <- as.factor(spreadsheet.data$Collection_end)
spreadsheet.data$Published <- as.factor(spreadsheet.data$Published)
spreadsheet.data$TraitID <- as.factor(spreadsheet.data$TraitID)

spreadsheet.data$Sex <- as.factor(spreadsheet.data$Sex)
spreadsheet.data$TraitID <- as.factor(spreadsheet.data$TraitID)
spreadsheet.data$Slope <- as.factor(spreadsheet.data$Slope)
spreadsheet.data$Drainage <- as.factor(spreadsheet.data$Drainage)
spreadsheet.data$Kingsolver_traits <- as.factor(spreadsheet.data$Kingsolver_traits)

## QUESTION 1 ECOLOGY

# First we make an "across" data frame for both slopes
## combine R2 and spreadsheet data
R2.data.among$TraitID <- as.factor(R2.data.among$TraitID)
data.for.models <- left_join(spreadsheet.data, R2.data.among,  by = "TraitID")

## filter by sex (because duplicates in 'Both')
data.for.models <- filter(data.for.models, Sex %in% c("M", "F"))

hist(data.for.models$R.2)  

# rename across slope set - because data.for.models will be used below too
data.for.models.across <- data.for.models

## now we make a South only data  frame
## combine R2 and spreadsheet data
R2.data.south$TraitID <- as.factor(R2.data.south$TraitID)
data.for.models.south <- left_join(spreadsheet.data, R2.data.south,  by = "TraitID")

## again, filter by sex to remove duplicates
data.for.models.south <- filter(data.for.models.south, Sex %in% c("M", "F"))

## here, filter by slope (doesn't really matter? Because ANOVA calculated only for South and TraitID not slope specific)
data.for.models.south <- data.for.models.south %>% filter(Slope == "South")

## remove duplicates so that we only have one R2 value per trait
data.for.models.south <- data.for.models.south[!duplicated(data.for.models.south$TraitID),]

## semi-join (return all values in across where there are matches in south)
## (this ensures that the TraitIDs are the same in both data frames)
data.for.models.across <- semi_join(data.for.models.across, data.for.models.south, by = "TraitID")

## remove duplicates in the across model, again so only 1 R2 per traitID
data.for.models.across <- data.for.models.across[!duplicated(data.for.models.across$TraitID),]

## QUESTION 2 EVOUTIONARY HISTORY

## first, make a data frame for caroni only
R2.data.caroni$TraitID <- as.factor(R2.data.caroni$TraitID)
data.for.models.caroni <- left_join(spreadsheet.data, R2.data.caroni,  by = "TraitID")

## filter by sex to remove duplicates from ('Both')
data.for.models.caroni <- filter(data.for.models.caroni, Sex %in% c("M", "F"))

## caroni only (again, idk if this really matters cuz ANOVA was run w only caroni and TraitID not specific 2 drainage)
data.for.models.caroni <- data.for.models.caroni %>% filter(Drainage == "Caroni")

## remove duplicates so only 1 R2 per traitID
data.for.models.caroni <- data.for.models.caroni[!duplicated(data.for.models.caroni$TraitID),]

## now we make a data frame for both caroni/oropuche
## combine R2 and spreadsheet data
R2.data.among.drainage$TraitID <- as.factor(R2.data.among.drainage$TraitID)
data.for.models.among.drainage <- left_join(spreadsheet.data, R2.data.among.drainage,  by = "TraitID")

## semi_join to return all values from both drainages that match values for caroni
## to make sure that the traitIDs are the same
data.for.models.among.drainage <- semi_join(data.for.models.among.drainage, data.for.models.caroni, by = "TraitID")

## remove duplicates so only 1 R2 per traitID
data.for.models.among.drainage <- data.for.models.among.drainage[!duplicated(data.for.models.among.drainage$TraitID),]

## QUESTION 3 INTRODUCTIONS

## first we will make intro only spreadsheet
## combine spreadsheet data w R2 for intro only ANOVA
R2.data.intro$TraitID <- as.factor(R2.data.intro$TraitID)
data.for.models.intro <- left_join(spreadsheet.data, R2.data.intro,  by = "TraitID")

## filter for M/F so no duplicates for 'Both'
data.for.models.intro <- data.for.models.intro %>% filter(Sex %in% c("M", "F"))

## only intro
data.for.models.intro <- filter(data.for.models.intro, Poptype == "Introduction")

## remove duplciates so only 1 R2 per traitID
data.for.models.intro <- data.for.models.intro[!duplicated(data.for.models.intro$TraitID),]

## now we make intro and natural
## same as above so can just use data.for.models
data.for.models.nat.intro <- data.for.models

## semi_join to return all values for both nat/intro that match traitID for intro
data.for.models.nat.intro <- semi_join(data.for.models.nat.intro, data.for.models.intro, by = "TraitID")

## remove duplicates so only 1 R2 per traitID
data.for.models.nat.intro <- data.for.models.nat.intro[!duplicated(data.for.models.nat.intro$TraitID),]


##%######################################################%##
#                                                          #
####                    GOOD MODELS                     ####
#                                                          #
##%######################################################%##

## note that I have excluded 'Other' for all of the trait models
## the values are super high so it's always significant and there have been some 
## issues w convergence

## full models

## we can just use data.for.models.across because that is for both slopes
full.trait.model <- glmer(R.2 ~ Kingsolver_traits + (1|Sex) + (1|StudyID), family = binomial, 
                          data = data.for.models.across[data.for.models.across$Kingsolver_traits %in% 
                                                          c("Other_morphology", "Size", 
                                                            "Physiology", "Behaviour", 
                                                            "Other_life_history", "Colour"),])  
full.sex.model <- glmer(R.2 ~ Sex + (1| Kingsolver_traits) +(1|StudyID), family = binomial, data = data.for.models.across)

summary(full.trait.model)
summary(full.sex.model)

library(car)
Anova(full.trait.model, type = "III")

## these we run just to check (so we can say that even tho we had 2 models it would have been the same w one)
full.trait.intxs <- glmer(R.2 ~ Kingsolver_traits*Sex + (1|StudyID), family = binomial,
                          data = data.for.models.across[data.for.models.across$Kingsolver_traits %in% 
                                                          c("Other_morphology", "Size", 
                                                            "Physiology", "Behaviour", 
                                                            "Other_life_history", "Colour"),])
full.trait.no.intxs <- glmer(R.2 ~ Kingsolver_traits + Sex + (1|StudyID), family = binomial,
                             data = data.for.models.across[data.for.models.across$Kingsolver_traits %in% 
                                                             c("Other_morphology", "Size", 
                                                               "Physiology", "Behaviour", 
                                                               "Other_life_history","Colour"),])

Anova(full.trait.intxs, type = "III")  # justify removing the interaction 
Anova(full.trait.no.intxs, type = "III")


summary(full.trait.model)
summary(full.trait.no.intxs)

data.for.models %>% 
  filter(Kingsolver_traits != "Other") %>% 
  ggplot(aes(x = Kingsolver_traits, y = R.2, color = Kingsolver_traits)) +
  geom_boxplot() +
  theme_bw()

data.for.models %>% 
  filter(Sex %in% c("M", "F")) %>% 
  ggplot(aes(x = Sex, y = R.2, color = Sex)) +
  geom_boxplot() +
  theme_bw()

## ecology question

## this is the same as above I did it twice so they could have diff names lol

hist(data.for.models.across$R.2)

across.model <- glmer(R.2 ~ Kingsolver_traits + Sex + (1|StudyID), family = binomial, 
                      control = glmerControl(optimizer = "bobyqa"),
                      data = data.for.models.across[data.for.models.across$Kingsolver_traits %in%
                                                      c("Other_morphology", "Size", 
                                                        "Physiology", "Behaviour", 
                                                        "Other_life_history", "Colour"),])

across.model.2 <- glmer(R.2 ~ Kingsolver_traits + (1|StudyID) + (1|Sex), family = binomial, 
                        control = glmerControl(optimizer = "bobyqa"),
                        data =data.for.models.across[data.for.models.across$Kingsolver_traits %in%
                                                       c("Other_morphology", "Size", 
                                                         "Physiology", "Behaviour", 
                                                         "Other_life_history", "Colour"),])

across.model.3 <- glmer(R.2 ~ Sex + (1|StudyID) + (1|Kingsolver_traits), family = binomial, 
                        control = glmerControl(optimizer = "bobyqa"),
                        data = data.for.models.across)

## rank deficient/doesn't converge
across.model.4 <- glmer(R.2 ~ Kingsolver_traits*Sex + (1|StudyID), family = binomial, 
                        data = data.for.models.across[data.for.models.across$Kingsolver_traits %in%
                                                        c("Other_morphology", "Size", 
                                                              "Physiology", "Behaviour", 
                                                              "Other_life_history", "Colour"),])

summary(across.model.2)
summary(across.model.3)

hist(data.for.models.south$R.2)


south.model <- glmer(R.2 ~ Kingsolver_traits + Sex + (1|StudyID), family = binomial, 
                     data = data.for.models.south[data.for.models.south$Kingsolver_traits %in% 
                                                    c( "Other_morphology", "Size", 
                                                       "Physiology", "Behaviour", 
                                                       "Other_life_history", "Colour"),])

south.model.2 <- glmer(R.2 ~ Kingsolver_traits + (1|StudyID) + (1|Sex), family = binomial, 
                       data =data.for.models.south[data.for.models.south$Kingsolver_traits %in% 
                                                     c( "Other_morphology", "Size", 
                                                        "Physiology", "Behaviour", 
                                                        "Other_life_history", "Colour"),])

south.model.3 <- glmer(R.2 ~ Sex + (1|StudyID) + (1|Kingsolver_traits), family = binomial, 
                       data = data.for.models.south)

soouth.model.4 <- glmer(R.2 ~ Kingsolver_traits*Sex + (1|StudyID), family = binomial, 
                        data = data.for.models.south[data.for.models.south$Kingsolver_traits %in% 
                                                       c( "Other_morphology", "Size", 
                                                          "Physiology", "Behaviour", 
                                                          "Other_life_history", "Colour"),])

summary(south.model.2)  # note the million warning messages???
summary(south.model.3)

data.for.models.south %>% 
  filter(Kingsolver_traits !="Other") %>% 
  ggplot(aes(x = Kingsolver_traits, y = R.2)) +
  geom_boxplot(aes(color = Kingsolver_traits)) +
  theme(legend.position = "none") +
  theme_bw()

data.for.models.south %>% 
  filter(Sex %in% c("M", "F")) %>% 
  ggplot(aes(x = Sex, y = R.2)) +
  geom_boxplot(aes(color = Sex)) +
  theme(legend.position = "none") +
  theme_bw()

ggplot(data = data.for.models.across[data.for.models.across$Kingsolver_traits %in% 
                                       c("Colour", "Behaviour", 
                                         "Other_morphology", "Other_life_history", 
                                         "Physiology", "Size"),],
       aes(y = R.2, color = "Across")) +
  geom_boxplot() +
  geom_boxplot(data = data.for.models.south[data.for.models.south$Kingsolver_traits %in% 
                                              c("Colour", "Behaviour", 
                                                "Other_morphology", "Other_life_history", 
                                                "Physiology", "Size"),]
               , aes(x = 1, y = R.2, color = "South"), position = position_dodge2()) +
  facet_wrap(~Kingsolver_traits) +
  theme_classic()

ggplot(data = data.for.models.across, aes(y = R.2, color = "Across")) +
  geom_boxplot() +
  geom_boxplot(data = data.for.models.south, 
               aes(x = 1, y = R.2, color = "South"), position = position_dodge2()) +
  theme_classic() +
  facet_wrap(~Sex)

## intro question

intro.model <- glmer(R.2 ~ Kingsolver_traits + Sex + (1|StudyID), family = binomial, 
                     data = data.for.models.intro[data.for.models.intro$Kingsolver_traits %in%
                                                    c("Other_morphology", "Size", 
                                                      "Physiology", "Behaviour", 
                                                      "Other_life_history", "Colour"),])

intro.model.2 <- glmer(R.2 ~ Kingsolver_traits + (1|StudyID) + (1|Sex), 
                       family = binomial, 
                       data =data.for.models.intro[data.for.models.intro$Kingsolver_traits %in% 
                                                     c("Other_morphology", "Size", 
                                                       "Physiology", "Behaviour", 
                                                       "Other_life_history", "Colour"),])

intro.model.3 <- glmer(R.2 ~ Sex + (1|StudyID) + (1|Kingsolver_traits), family = binomial, 
                       data = data.for.models.intro)

intro.model.4 <- glmer(R.2 ~ Sex*Kingsolver_traits + (1|StudyID), family = binomial, 
                       data = data.for.models.intro[data.for.models.intro$Kingsolver_traits %in%
                                                      c("Other_morphology", "Size", 
                                                        "Physiology", "Behaviour", 
                                                        "Other_life_history", "Colour"),])

summary(intro.model.2)
summary(intro.model.3)

ggplot(data = data.for.models.across, aes(y = R.2, color = "Both")) +
  geom_boxplot() + 
  geom_boxplot(data = data.for.models.intro, aes(x = 1, y = R.2, color = "Intro")) +
  facet_wrap(~Kingsolver_traits) + theme_classic()

ggplot(data = data.for.models.across, aes(x = 0, y = R.2, color = "Both")) +
  geom_boxplot() + 
  geom_boxplot(data = data.for.models.intro, aes(x = 1, y = R.2, color = "Intro")) +
  theme_classic() +
  facet_wrap(~Sex)

## evolutionary history question

## there are some issues w convergence in caroni.model.2 but these go away if we exclude colour

among.drainage.model <- glmer(R.2 ~ Kingsolver_traits + Sex + (1|StudyID), family = binomial, 
                              data = data.for.models.among.drainage[data.for.models.among.drainage$Kingsolver_traits %in% 
                                                                      c( "Other_morphology", "Size", 
                                                                         "Physiology", "Behaviour", 
                                                                         "Other_life_history", "Colour"),])

among.drainage.model.2 <- glmer(R.2 ~ Kingsolver_traits + (1|StudyID) + (1|Sex), 
                                family = binomial, 
                                data =data.for.models.among.drainage[data.for.models.among.drainage$Kingsolver_traits %in% 
                                                                       c( "Other_morphology", "Size", 
                                                                          "Physiology", "Behaviour", 
                                                                          "Other_life_history", "Colour"),])

among.drainage.model.3 <- glmer(R.2 ~ Sex + (1|StudyID) + (1|Kingsolver_traits), family = binomial, data = data.for.models.among.drainage)

among.drainage.model.4 <- glmer(R.2 ~ Kingsolver_traits*Sex + (1|StudyID), family = binomial, 
                              data = data.for.models.among.drainage[data.for.models.among.drainage$Kingsolver_traits %in% 
                                                                      c("Other_morphology", "Size", 
                                                                         "Physiology", "Behaviour", 
                                                                         "Other_life_history", "Colour"),])

summary(among.drainage.model.2)
summary(among.drainage.model.3)

caroni.model <- glmer(R.2 ~ Kingsolver_traits + Sex + (1|StudyID), family = binomial, 
                      data = data.for.models.caroni[data.for.models.caroni$Kingsolver_traits %in% 
                                                      c("Other_morphology", "Size",
                                                        "Physiology",  "Behaviour", 
                                                        "Other_life_history",  "Colour"),])

caroni.model.2 <- glmer(R.2 ~ Kingsolver_traits + (1|StudyID) + (1|Sex), 
                        family = binomial, 
                        data =data.for.models.caroni[data.for.models.caroni$Kingsolver_traits %in% 
                                                       c("Other_morphology", "Size",
                                                         "Physiology",  "Behaviour", 
                                                         "Other_life_history",  "Colour"),])

caroni.model.3 <- glmer(R.2 ~ Sex + (1|StudyID) + (1|Kingsolver_traits), family = binomial, data = data.for.models.caroni)

caroni.model.4 <- glmer(R.2 ~ Kingsolver_traits*Sex + (1|StudyID), family = binomial, 
                      data = data.for.models.caroni[data.for.models.caroni$Kingsolver_traits %in% 
                                                      c("Other_morphology", "Size",
                                                        "Physiology",  "Behaviour", 
                                                        "Other_life_history",  "Colour"),])

summary(caroni.model.2)
summary(caroni.model.3)

data.for.models.among.drainage %>% 
  filter(Kingsolver_traits != "Other") %>% 
  ggplot(aes(x = Kingsolver_traits, y = R.2)) +
  geom_boxplot(aes(color = Kingsolver_traits)) +
  geom_jitter(size = 4, alpha = 0.1, aes(color = Kingsolver_traits)) +
  theme(legend.position = "none") +
  theme_classic()

data.for.models.caroni %>% 
  filter(Kingsolver_traits != "Other") %>% 
  ggplot(aes(x = Kingsolver_traits, y = R.2, color = Kingsolver_traits)) +
  geom_boxplot() +
  geom_jitter(size = 4, alpha = 0.1) +
  theme(legend.position = "none") +
  theme_classic()

ggplot(data = data.for.models.among.drainage, aes(y = R.2, color = "Both drainages")) +
  geom_boxplot() + 
  geom_boxplot(data = data.for.models.caroni, aes(x = 1, y = R.2, color = "Caroni only")) +
  theme_classic() +
  facet_wrap(~Sex)

ggplot(data = data.for.models.among.drainage, aes(y = R.2, color = "Both drainages")) +
  geom_boxplot() + 
  geom_boxplot(data = data.for.models.caroni, aes(x = 1, y = R.2, color = "Caroni only")) +
  facet_grid(ncol = vars(Kingsolver_traits), nrow = vars(Sex)) +
  theme_classic()


#### here down is the old stuff ####

# set working directory
# setwd("")
# wd<- getwd()

# Import and tidy
spreadsheet.data <- read.csv(paste(wd,'/Data/MetaData.csv',sep=""), header=TRUE, sep=",")
R2.data <- read.csv(paste(wd,'/Data/TraitR2.csv',sep=""), header=TRUE, sep=",")

# spreadsheet.data is the data extracted for the meta-analysis
# R2.data is the output of the ANOVA loop


# AH 2020-10-27 deleted because no longer diff between spreadsheets
# names(R2.data)[names(R2.data) == "TraitID"] <- "sex_TraitID"  # so same in both spreadsheets
# all other sex_TraitIDs below changed as well

str(spreadsheet.data)
str(R2.data)

spreadsheet.data$StudyID <- as.factor(spreadsheet.data$StudyID)
spreadsheet.data$Collection_start <- as.factor(spreadsheet.data$Collection_start)
spreadsheet.data$Collection_end <- as.factor(spreadsheet.data$Collection_end)
spreadsheet.data$Published <- as.factor(spreadsheet.data$Published)
spreadsheet.data$TraitID <- as.factor(spreadsheet.data$TraitID)

R2.data$TraitID <- as.factor(R2.data$TraitID)

# This (data.for.models) is the data to use
data.for.models <- left_join(spreadsheet.data, R2.data,  by = "TraitID")
hist(data.for.models$R.2)  

data.for.models$Sex <- as.factor(data.for.models$Sex)
data.for.models$TraitID <- as.factor(data.for.models$TraitID)
data.for.models$Slope <- as.factor(data.for.models$Slope)
data.for.models$Drainage <- as.factor(data.for.models$Drainage)

str(data.for.models)

# Collapsing traits
data.for.models$TraitType2[data.for.models$TraitType2 == "Diet"] <- "Other"
data.for.models$TraitType2[data.for.models$TraitType2 == "Physiology"] <- "Other"
data.for.models$TraitType2[data.for.models$TraitType2 == "Life history"] <- "Other"

hist(data.for.models$R.2)

# data.for.models <- data.for.models[!duplicated(data.for.models$TraitID),]
# this line gets rid of duplicates - but also lots of detail

##%######################################################%##
#                                                          #
####                  EXPLORING DATA                    ####
#                                                          #
##%######################################################%##


plot(table(data.for.models$R.2), ylab = "Frequency", xlab = "R.2",
     main = "Frequency of R2 values", cex.lab = 1.5, cex.main = 1.5)

data.for.models %>% filter(StudyType == "Wildcaught" & Sex == "Both" & Drainage != "NA") %>% 
  ggplot(aes(x = TraitType2, y = R.2)) +
  geom_boxplot() +
  facet_wrap(~Drainage)

data.for.models %>% filter(StudyType == "Wildcaught" & Sex == "Both" & Slope != "Other") %>% 
  ggplot(aes(x = TraitType2, y = R.2)) +
  geom_boxplot() +
  facet_wrap(~Slope)

(cgtally <-
  data.for.models %>% 
  filter(StudyType %in% c("Common Garden (F1)", ("Common Garden (F2)")) & Slope %in% c("North", "South")) %>% 
  group_by(Slope, Sex, TraitType2) %>% tally())

(wctally <-
  data.for.models %>% 
  filter(StudyType == "Wildcaught" & Slope %in% c("North", "South") & Sex %in% c("M", "F", "Both")) %>% 
  group_by(Slope, Sex, TraitType2) %>% tally())

(yeartally <-
    data.for.models %>% 
    filter(StudyType == "Wildcaught" & Slope %in% c("North", "South") & Sex %in% c("M", "F", "Both")) %>% 
    group_by(Collection_end, StudyID) %>% tally())

## MOD 1 Traits/both sexes

# WC 
(overall.trait.plot.WC <-  # many more studies in the south....  
    data.for.models %>% filter(StudyType == "Wildcaught" & Sex == "Both") %>% 
    ggplot(aes(x = TraitType2, y = R.2)) +
    theme_bw() + 
    geom_violin(aes(fill = TraitType2)) + 
    geom_jitter(width = 0.2, alpha = 0.3) +
    labs(x = "Trait type", y = "R2", title = "Both sexes, Wildcaught"))

# CG
(overall.trait.plot.WC <-  # only studies in the south  
    data.for.models %>% filter(StudyType 
                               %in% c("Common Garden (F1)", "Common Garden (F2)") 
                               & Sex == "Both") %>% 
    ggplot(aes(x = TraitType2, y = R.2)) +
    theme_bw() + 
    geom_violin(aes(fill = TraitType2)) + 
    labs(x = "Trait type", y = "R2", title = "Both sexes, Common Garden"))

## MOD 2 Sex-specific

# WC

(male.trait.plot.WC <- 
data.for.models %>% filter(StudyType == "Wildcaught" & Sex == "M") %>% 
  ggplot(aes(x = TraitType2, y = R.2)) +
  theme_bw() + 
    geom_boxplot(aes(fill = Slope)) + 
  labs(x = "Trait type", y = "R2", title = "Males WC"))

# CG
(male.trait.plot.CG <- 
  data.for.models %>% filter(StudyType 
                             %in% c("Common Garden (F1)", "Common Garden (F2)") 
                             & Sex == "M") %>% 
  ggplot(aes(x = TraitType2, y = R.2)) +
  theme_bw() +
    geom_violin(aes(fill = TraitType2)) +
  labs(x = "Trait type", y = "R2", title = "Males CG"))

## 

# WC

(female.trait.plot.WC <-
    data.for.models %>% filter(StudyType == "Wildcaught" & Sex == "F") %>% 
    ggplot(aes(x = TraitType2, y = R.2)) +
    theme_bw() +
    geom_violin(aes(fill = TraitType2)) +
    labs(x = "Trait type", y = "R2", title = "Females WC"))

# CG
(female.trait.plot.CG <-
  data.for.models %>% filter(StudyType 
                             %in% c("Common Garden (F1)", "Common Garden (F2)") 
                             & Sex == "F") %>% 
  ggplot(aes(x = TraitType2, y = R.2)) +
  theme_bw() +
    geom_violin(aes(fill = TraitType2)) +
  labs(x = "Trait type", y = "R2", title = "Females CG"))

## MOD 4 Collection_End
data.for.models %>% filter(StudyType == "Wildcaught" 
                           & Sex == "Both" 
                           & Collection_end != "NA") %>% 
  ggplot(aes(x = Collection_end, y = R.2)) +
  theme_bw() +
  geom_violin(aes(fill = Collection_end)) +
  labs(x = "Collection end year", y = "R2", title = "Last year of data collection")

data.for.models %>% filter(StudyType == "Wildcaught" 
                           & Sex == "Both" 
                           & Collection_end != "NA") %>% 
  ggplot(aes(x = Collection_end, y = R.2)) +
  theme_bw() +
  stat_summary(fun.y = mean, color = "deepskyblue2") +
  labs(x = "Collection end year", y = "R2", title = "Average R2 per collection end year") +
  ylim(c(0,1))

# Correlation plots! 

(tbl1 <- table(data.for.models$TraitType2, data.for.models$Slope))
chisq.test(tbl1)  

(tbl2 <- table(data.for.models$TraitType2, data.for.models$Paired))
chisq.test(tbl2)  

(tbl3 <- table(data.for.models$Paired, data.for.models$Slope))
chisq.test(tbl3)  

corrplot(tbl1, is.cor = FALSE)
corrplot(tbl2, is.cor = FALSE)
corrplot(tbl3, is.cor = FALSE)

(tbl4 <- table(data.for.models$TraitType2, data.for.models$Sex))
tbl4
corrplot(tbl4, is.cor = FALSE)

## These below are for regression to the mean

data.for.models %>% 
  filter(Collection_end != "NA") %>% 
  group_by(StudyID, Collection_end) %>% 
  ggplot(aes(x = R.2)) + geom_histogram(binwidth = 0.03) 

data.reg2m <- data.for.models %>% 
  filter(Collection_end !="NA") %>% 
  group_by(Collection_end, StudyID) %>% 
  summarise(meanR2 = mean(R.2))

head(data.reg2m)

# individual R2 for Collection_end years 
# don't know what to do because not paired ... 

R2 <- data.reg2m$meanR2
names(R2) <- paste(data.reg2m$Collection_end)
print(R2)
mean(R2, na.rm = TRUE)  

with(data.reg2m, plot(x = Collection_end, y = R2))

(ggplot(data.reg2m, aes(x = meanR2)) + 
    geom_density() +
    geom_vline(aes(xintercept = mean(meanR2)), 
               color = "blue",
               linetype = "dashed",
               size = 1))

ggplot(data.reg2m, aes(x = Collection_end, y = meanR2)) +
  geom_point() +
  theme_classic() +
  stat_summary(aes(group = 1), fun = mean, colour = "red", geom = "line") +
  geom_hline(yintercept = 0.3467688, colour = "blue", linetype = "dashed")

data.for.models %>%
  ggplot(aes(x = R.2, group = TraitType2, color = TraitType2)) +
  geom_density(alpha = 0.6) + theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))


##%######################################################%##
#                                                          #
####                   DISTRIBUTION                     ####
#                                                          #
##%######################################################%##

#log can't have 0s, so diff dataset
dist.data.log<-R2.data$R.2+.001
dist.data<-R2.data$R.2
plot(dist.data, pch=20)
#empirical data
plotdist(dist.data, histo = TRUE, demp = TRUE)

descdist(dist.data, discrete=FALSE, boot=500)
descdist(dist.data.log, discrete=FALSE, boot=500)

fit_b<-fitdist(dist.data, "beta", method="mme")
denscomp(fit_b, ylim=(1))
cdfcomp (fit_b,ylim=(1))
qqcomp  (fit_b)
ppcomp  (fit_b)


fit_bn<-fitdist(dist.data, "nbinom", method="mme")
denscomp(fit_bn,  ylim=(1))
cdfcomp (fit_bn,  ylim=(1))
qqcomp  (fit_bn ) #I cannot plot these two and I have no idea why
ppcomp  (fit_bn)

#diff datasets can't be compared to each other
fit_ln<-fitdist(dist.data.log, "lnorm", method="mme")
cdfcomp(fit_ln, xlogscale = TRUE, ylogscale = TRUE) 
denscomp(fit_ln,  ylim=(1))
qqcomp  (fit_ln )
ppcomp  (fit_ln)
#idk if I'm doing something wrong or this is a bad fit

gofstat(fit_b)
gofstat(fit_bn)
gofstat(fit_ln)

par(mfrow=c(2,2))
plot.legend <- c("beta", "nbinom")
denscomp(list(fit_b, fit_bn), legendtext = plot.legend, ylim=(1))
cdfcomp (list(fit_b, fit_bn), legendtext = plot.legend, ylim=(1))
qqcomp  (list(fit_b, fit_bn), legendtext = plot.legend)
ppcomp  (list(fit_b, fit_bn), legendtext = plot.legend)
par(mfrow = c(1,1))

##%######################################################%##
#                                                          #
####                    GOOD MODELS                     ####
#                                                          #
##%######################################################%##

# I think that the best model overall is sex.mod6!! #

### QUESTION 1. IS THERE A DIFFERENCE OVERALL BETWEEN HIGH/LOW ###

# I think this one works THE BEST for Q1 -- BUT only South
mod1.south <- glmer(R.2 ~ TraitType2 + (1|StudyID), 
                       family = binomial,
              data = data.for.models[
                data.for.models$StudyType == "Wildcaught" 
                & data.for.models$Sex == "Both"
                & data.for.models$Slope == "South",
                ])
summary(mod1.south)
AIC(mod1.south)  # 441.87

### QUESTION 2. IS PARALLELISM DIFFERENT BETWEEN THE SEXES? ###

# 2020-10-14 omitting females from this formula
# (there is the weird one where they added testonerone to see change in colour)
# (I think it is making colour significant when it has not been)
# (^ This is studyID 30, so trying to omit)


sex.mod6 <- glmer(R.2 ~ Sex + TraitType2 + Slope + (1|StudyID),
                  family = binomial, 
                  data =
                    data.for.models[data.for.models$StudyType == "Wildcaught"
                                    & data.for.models$Sex %in% c("M", "F"),])
summary(sex.mod6)
AIC(sex.mod6)  #2162.201 Oct 28th (AP)

## with slope type column (this indicates if purely south or purely north or mixed)

sex.mod7 <- glmer(R.2 ~ Sex + TraitType2 + Slope.type + (1|StudyID),
                  family = binomial, 
                  data =
                    data.for.models[data.for.models$StudyType == "Wildcaught"
                                    & data.for.models$Sex %in% c("M", "F"),])
summary(sex.mod7)
AIC(sex.mod7)  #2152.358 Dec 2nd (AP)
tab_model(sex.mod7)


#### QUESTION 3 - IS THERE A DIFFERENCE BETWEEN THE SLOPES? ####
# For this question, using only "Both" sexes #

slope.mod1 <- glmer(R.2 ~ Slope + (1|StudyID),
                     family = binomial, 
                     data = data.for.models[data.for.models$StudyType == "Wildcaught"
                                            & data.for.models$Sex == "Both",])
summary(slope.mod1)
AIC(slope.mod1)  # 695.6815

#### QUESTION 4. IS THERE REGRESSION TOWARDS THE MEAN ####
## Here, we are looking at both sexes and Collection_end ##

# updated so just a glm! 2020-09-29


# As per the advice of Guillaume from QCBS - switching to numeric.
data.for.models$Collection_end <- as.numeric(paste(data.for.models$Collection_end))
str(data.for.models)

time.mod2 <- glm(R.2 ~ Collection_end,
                   family = binomial,
                   data = data.for.models[data.for.models$Sex %in% c("M", "F") 
                                          & data.for.models$StudyType == "Wildcaught",])
summary(time.mod2)
confint(time.mod2)


# this one runs, but probably shouldn't ignore Study.ID?
time.mod3 <- glm(R.2 ~ Collection_end, 
                 family = binomial, 
                 data = data.for.models[data.for.models$Sex == "Both" 
                                        & data.for.models$StudyType == "Wildcaught",])
summary(time.mod3)


##%######################################################%##
#                                                          #
####                   OTHER MODELS                     ####
#                                                          #
##%######################################################%##

## Q1 

# has drainages, won't converge - crap
mod1.drain <- glmer(R.2 ~ TraitType2 + (1|StudyID/Drainage), 
                    family = binomial, 
                    data = data.for.models[
                      data.for.models$StudyType == "Wildcaught" 
                      & data.for.models$Sex == "Both",])

summary(mod1.drain)

# Again, this one won't converge - crappy crap
mod1.slope <- glmer(R.2 ~ TraitType2 + (1|StudyID/Slope), 
                    family = binomial, 
                    data = data.for.models[
                      data.for.models$StudyType == "Wildcaught" 
                      & data.for.models$Sex == "Both",])

summary(mod1.slope)
AIC(mod1.slope)



# Only North - Won't run (error) 
mod1.north <- glmer(R.2 ~ TraitType2 + (1|StudyID), 
                    family = binomial,
                    data = data.for.models[
                      data.for.models$StudyType == "Wildcaught" 
                      & data.for.models$Sex == "Both"
                      & data.for.models$Slope == "North",
                      ])
summary(mod1.north)

# Also doesn't run
mod1b.north <- glmer(R.2 ~ TraitType2 + (1|Study.ID), 
                     family = binomial,
                     data = data.for.models[
                       data.for.models$StudyType == "Wildcaught" 
                       & data.for.models$Sex == "Both"
                       & data.for.models$Slope == "North",
                       ])
summary(mod1b.north)

## Q2

sex.mod1 <- glmer(R.2 ~ Sex*Slope + (1|StudyID),
                  family = binomial, 
                  data = data.for.models[data.for.models$StudyType == "Wildcaught"
                                         & data.for.models$Sex %in% c("M", "F"),])
summary(sex.mod1)
AIC(sex.mod1) # 1995.517

sex.mod3 <- glmer(R.2 ~ Sex + (1|Study.ID),
                  family = binomial, 
                  data = data.for.models[data.for.models$StudyType == "Wildcaught"
                                         & data.for.models$Sex %in% c("M", "F"),])
summary(sex.mod3)
AIC(sex.mod3)  # 2387.036

# First, I am doing both "M" and "F" (i.e. not "both"Both", and not separate)

# This model probably makes the most sense out of the broad sex specific ones 
# sex*slope also runs (w other models down below)
sex.mod2 <- glmer(R.2 ~ Sex+Slope + (1|StudyID),
                  family = binomial, 
                  data = data.for.models[data.for.models$StudyType == "Wildcaught"
                                         & data.for.models$Sex %in% c("M", "F"),])
summary(sex.mod2)
AIC(sex.mod2)  # 1994.825

# SECOND - specific looking at sex and traits: 
# Again, here both "M" and "F"
# In this section I am NOT including colour or life history
# Colour = only males
# Life history = only females

# This one??
sex.mod4 <- glmer(R.2 ~ Sex*TraitType2 + (1|StudyID),
                  family = binomial, 
                  data =
                    data.for.models[data.for.models$StudyType == "Wildcaught"
                                    & data.for.models$Sex %in% c("M", "F"),])
summary(sex.mod4)
AIC(sex.mod4)

sex.mod5 <- glmer(R.2 ~ Sex + TraitType2 + (1|StudyID),
                  family = binomial, 
                  data =
                    data.for.models[data.for.models$StudyType == "Wildcaught"
                                    & data.for.models$Sex %in% c("M", "F"),])
summary(sex.mod5)
AIC(sex.mod5)


# Singular fit
sex.mod7.wc <- glmer(R.2 ~ Sex + (1|Study.ID/Slope),
                     family = binomial, 
                     data =
                       data.for.models[data.for.models$StudyType == "Wildcaught"
                                       & data.for.models$Sex %in% c("M", "F")
                                       & data.for.models$TraitType2 %in% c("Diet", "Morphometric", "Other", "Physiology", "Behaviour"),])
summary(sex.mod7.wc)
AIC(sex.mod7.wc) 

# Now, looking at males and females separately (i.e. only "M" or only "F")

# Males only ... 

trait.mod3.males <- glmer(R.2 ~ Slope + TraitType2 + (1|Study.ID), 
                          family = binomial,
                          data = data.for.models[
                            data.for.models$StudyType == "Wildcaught" 
                            & data.for.models$Sex == "M",
                            ])

summary(trait.mod3.males)
AIC(trait.mod3.males)  # 1292.874

trait.mod5.males <- glmer(R.2 ~ TraitType2 + (1|Study.ID), 
                          family = binomial,
                          data = data.for.models[
                            data.for.models$StudyType == "Wildcaught" 
                            & data.for.models$Sex == "M",
                            ])

summary(trait.mod5.males)
AIC(trait.mod5.males)  # 1495.068

# The next model is male only, and ONLY SOUTH

trait.mod6.males <- glmer(R.2 ~ TraitType2 + (1|Study.ID), 
                          family = binomial,
                          data = data.for.models[
                            data.for.models$StudyType == "Wildcaught" 
                            & data.for.models$Sex == "M"
                            & data.for.models$Slope == "South",
                            ])

summary(trait.mod6.males)
AIC(trait.mod6.males)  # 806.9988

# wont converge
trait.mod7.males <- glmer(R.2 ~ TraitType2 + (1|Study.ID), 
                          family = binomial,
                          data = data.for.models[
                            data.for.models$StudyType == "Wildcaught" 
                            & data.for.models$Sex == "M"
                            & data.for.models$Slope == "North",
                            ])

summary(trait.mod7.males)

# Q4 with "Both"

time.mod1 <- glm(R.2 ~ Collection_end,
                 family = binomial,
                 data = data.for.models[data.for.models$Sex == "Both" 
                                        & data.for.models$StudyType == "Wildcaught",])
summary(time.mod1)

##%######################################################%##
#                                                          #
####                   ALEXIS' FIGS                     ####
#                                                          #
##%######################################################%##

# these are for my term paper

# 1 ecology (slope)
data.for.models %>% filter(StudyType == "Wildcaught" & Sex %in% c("M", "F") & Slope %in% c("North", "South")) %>% 
  ggplot(aes(x = TraitType2, y = R.2)) +
  geom_boxplot() +
  facet_wrap(~Slope)

# 2 evolutionary history
data.for.models %>% filter(StudyType == "Wildcaught" & Sex %in% c("M", "F") & Drainage %in% c("Caroni", "Oropuche")) %>% 
  ggplot(aes(x = TraitType2, y = R.2, color = TraitType2)) +
  geom_violin(size = 1.25) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text = element_text(size = 12)) +
  facet_wrap(~Drainage, ncol = 1)

# 3 time [NA]
# 4 trait types

data.for.models %>% filter(Sex %in% c("M", "F") & Slope %in% c("North", "South")) %>% 
  ggplot(aes(x = R.2, y = TraitType2, color = TraitType2)) +
  geom_boxplot(size = 1.25) + 
  geom_jitter(width = 0.1, alpha = 0.2) +
  xlab("R2") + ylab("Trait type \n") +
  xlim(0,1) +
  theme(legend.position = "none") +
  theme_classic() +
  facet_wrap(~Sex)

data.for.models %>% filter(Sex %in% c("M", "F") & Slope %in% c("North", "South")) %>% 
  ggplot(aes(x = R.2, y = TraitType2, color = TraitType2)) +
  geom_boxplot(size = 1.25) + 
  geom_jitter(width = 0.1, alpha = 0.2) +
  xlab("R2") + ylab("Trait type \n") +
  xlim(0,1) +
  theme(legend.position = "none") +
  theme_classic() +
  facet_wrap(~Slope)

data.for.models %>%
  ggplot(aes(R.2, color = TraitType2)) +
  geom_freqpoly() + theme_classic() + 
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text  = element_text(size = 12),
        legend.title = element_text(size = 14)) 

# 5 sex
data.for.models %>% filter(StudyType == "Wildcaught" & Sex %in% c("M", "F")) %>% 
  ggplot(aes(x = Sex, y = R.2, color = Sex)) +
  geom_boxplot(size = 1.25) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

# 6. wc v cg

data.for.models %>% filter(Sex %in% c("M", "F")) %>% 
  ggplot(aes(x = StudyType, y = R.2)) + geom_boxplot() +
