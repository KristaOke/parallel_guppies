# GLMMs for parallel_guppies! 
# 2021-07-26 - idk what happened but this was replaced by an old version
# going to copy/paste from working script (new)
# pushed old stuff to bottom

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
library(olsrr)
library(ggpubr)
library(PupillometryR)
library(MuMIn)


# set working directory
# setwd("")
# wd<- getwd()


# import and tidy 

## These are all the updated sheets on the Drive
## spreadsheet.data is the data extracted for the meta-analysis
## R2.data.among are the output of the ANOVA loops

spreadsheet.data <- read.csv(paste(wd,'/Data/MetaData.csv',sep=""), header=TRUE, sep=",")
R2.data.among <- read.csv(paste(wd,'/Data/TraitR2_among.csv',sep=""), header=TRUE, sep=",")
R2.data.south <- read.csv(paste(wd,'/Data/TraitR2_south.csv',sep=""), header=TRUE, sep=",")
R2.data.intro <- read.csv(paste(wd,'/Data/TraitR2_intro.csv',sep=""), header=TRUE, sep=",")
R2.data.caroni <- read.csv(paste(wd,'/Data/TraitR2_Caroni.csv',sep=""), header=TRUE, sep=",")
R2.data.among.drainage <- read.csv(paste(wd,'/Data/TraitR2_Among_Drainage.csv',sep=""), header=TRUE, sep=",")
R2.data.intro.broad <- read.csv(paste(wd, '/Data/TraitR2_intro_broad.csv', sep = ""), header = TRUE, sep = "")

## fix structure
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

R2.data.among$TraitID <- as.factor(R2.data.among$TraitID)
R2.data.south$TraitID <- as.factor(R2.data.south$TraitID)
R2.data.intro$TraitID <- as.factor(R2.data.intro$TraitID)
R2.data.intro.broad$TraitID <- as.factor(R2.data.intro.broad$TraitID)
R2.data.caroni$TraitID <- as.factor(R2.data.caroni$TraitID)
R2.data.among.drainage$TraitID <- as.factor(R2.data.among.drainage$TraitID)

# combine R2 and spreadsheet data
## prep for when we bind them together, this variable will go in the model to indicate the type of R2
R2.data.among$method <- "all"
R2.data.south$method <- "south"
R2.data.intro$method <- "only_natural"
R2.data.intro.broad$method <- "only_natural_broad"
R2.data.caroni$method <- "caroni"
R2.data.among.drainage$method <- "both.drainages"

## get relevant data with one entry for each Trait
### I was inclusive with columns, many probably aren't needed so you can cut them out
data.all <- inner_join(spreadsheet.data, R2.data.among, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:52)%>% #I selected only the columns that have information that applies at the trait level
  distinct(TraitID, .keep_all = TRUE) #removes replicated TraitIDs and retains the columns

## repeat for south only R2 and others...
data.south <- inner_join(spreadsheet.data, R2.data.south, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:52)%>% 
  distinct(TraitID, .keep_all = TRUE) 

data.intro<- inner_join(spreadsheet.data, R2.data.intro, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:52)%>% 
  distinct(TraitID, .keep_all = TRUE) 

data.intro.broad<- inner_join(spreadsheet.data, R2.data.intro.broad, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:52)%>% 
  distinct(TraitID, .keep_all = TRUE) 

data.caroni<- inner_join(spreadsheet.data, R2.data.caroni, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:52)%>% 
  distinct(TraitID, .keep_all = TRUE) 

data.among.drainage<- inner_join(spreadsheet.data, R2.data.among.drainage, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:52)%>% 
  distinct(TraitID, .keep_all = TRUE)

## then we can bind them together as we wish! First ecology...
data.for.ecology.models<-rbind(data.all,data.south) %>% 
  arrange(TraitID) %>% #puts them in a nice order
  group_by(TraitID) %>% #groups them for the count
  filter(n() > 1) #filters only trait IDs that have more than 1 entry within the group (n = 204 traits) %>% 

ecology.data.males <- data.for.ecology.models %>% 
  filter(Sex == "M" & Kingsolver_traits !="Other") %>% 
  ungroup(TraitID)

ecology.data.females <- data.for.ecology.models %>% 
  filter(Sex == "F" & Kingsolver_traits !="Other") %>% 
  ungroup(TraitID)

## evolutionary history question
data.for.evolhist.models<-rbind(data.caroni,data.among.drainage) %>% 
  arrange(TraitID) %>% #puts them in a nice order
  group_by(TraitID) %>% #groups them for the count
  filter(n() > 1) #filters only trait IDs that have more than 1 entry within the group (n = 204 traits)

evolhist.data.males <- data.for.evolhist.models %>% 
  filter(Sex == "M" & Kingsolver_traits !="Other") %>% 
  ungroup(TraitID)
evolhist.data.females <- data.for.evolhist.models %>% 
  filter(Sex == "F" & Kingsolver_traits !="Other") %>% 
  ungroup(TraitID)

## Intro "strict", I'm pretty sure these are the right csvs for this question
data.for.intro.models<-rbind(data.all,data.intro) %>% 
  arrange(TraitID) %>% #puts them in a nice order
  group_by(TraitID) %>% #groups them for the count
  filter(n() > 1) #filters only trait IDs that have more than 1 entry within the group (n = 204 traits)

time.data.males <- data.for.intro.models %>% 
  filter(Sex == "M" & Kingsolver_traits !="Other") %>% 
  ungroup(TraitID)
time.data.females <- data.for.intro.models %>% 
  filter(Sex == "F" & Kingsolver_traits !="Other") %>% 
  ungroup(TraitID)


## Intro "broad", 
data.for.intro.models.broad<-rbind(data.all,data.intro.broad) %>% 
  arrange(TraitID) %>% #puts them in a nice order
  group_by(TraitID) %>% #groups them for the count
  filter(n() > 1) #filters only trait IDs that have more than 1 entry within the group (n = 204 traits)

### as a note and fail safe I would run the grouping with TraitID on the sex specific dfs to make sure theres two traitID entries for each

# Overall models ----

## remove 'both'
### (Because when calculated by hand, I would do M/F/Both with same data)
data.all <- data.all %>% filter(Sex %in% c("M", "F"))  

mean(data.all$R.2)
sd(data.all$R.2)
#sd(data.all.no.colour$R.2)

## remove other (because model below will not converge with other)
data.all.traits <- data.all %>% filter(!Kingsolver_traits == "Other")

## creating a dataset w no colour (to compare w those w colour)
data.all.no.colour <- data.all %>% filter(!Kingsolver_traits == "Colour")

## overall traits model ----

## All traits - this WILL NOT RUN with other included
(all.model.traits <- glmer(R.2 ~ Kingsolver_traits +  (1|StudyID), 
                           data = data.all.traits, family = binomial)) %>% summary()

## Rearing enviro model ----
data.all.rear <- data.all.traits %>% filter(StudyType %in% c("Common Garden (F2)", "Wildcaught"))  # won't run w CG F1 (not a lot anyway)
(all.model.rearing <- glmer(R.2 ~ StudyType +  (1|StudyID), data = data.all.rear, family = binomial)) %>% summary()

##  overall sex model ----
(all.model.sex <- glmer(R.2 ~ Sex + (1|StudyID), data = data.all, family = binomial)) %>% summary()

## Removed colour 
(all.model.sex.no.colour <- glmer(R.2 ~ Sex + (1|StudyID), data = data.all.no.colour, family = binomial)) %>% summary()

# multivariate models ----

(sex.and.triats <- glmer(R.2 ~ Kingsolver_traits + Sex + (1|StudyID), data = data.all, family = binomial)) %>% summary()
Anova(sex.and.triats, type = "II")

(sex.and.rear <- glmer(R.2 ~ StudyType + Sex + (1|StudyID), data = data.all.rear, family = binomial)) %>% summary()

Anova(sex.and.rear, type = 2)

# Traits descriptive stuff for in text ----

(Colourtrait <- data.all %>% filter(Kingsolver_traits == "Colour")) %>% summary()
(LHtrait <- data.all %>% filter(Kingsolver_traits == "Other_life_history")) %>% summary()
(Sizetrait <- data.all %>% filter(Kingsolver_traits == "Size")) %>% summary()
(Morphtrait <- data.all %>% filter(Kingsolver_traits == "Other_morphology")) %>% summary()
(Othertrait <- data.all %>% filter(Kingsolver_traits == "Other")) %>% summary()
(Phystrait <- data.all %>% filter(Kingsolver_traits == "Physiology")) %>% summary()

# Rearing enviro for in-text
(cg <- data.all %>% filter(StudyType == "Common Garden (F2)")) %>% summary()
(wc <- data.all %>% filter(StudyType == "Wildcaught")) %>% summary()

# Sex
(Malewithcolour <- data.all %>% filter(Sex == "M")) %>% summary()
(Malenocolour <- data.all.no.colour %>% filter(Sex == "M")) %>% summary()
(Femalesex <- data.all %>% filter(Sex == "F")) %>% summary()


# Determinants (ecological complex, evol hist, contemp evol) ----

## Sample size tables ----

### I don't think that you can just tally by method

## Ecology
with(data.for.ecology.models, table(Kingsolver_traits, Sex))
with(data.for.ecology.models, table(Kingsolver_traits))
with(data.for.ecology.models, table(Sex))
data.for.ecology.models %>% group_by(StudyID) %>% tally()

## Intro
with(data.for.intro.models, table(Kingsolver_traits, Sex))
with(data.for.intro.models, table(Kingsolver_traits))
with(data.for.intro.models, table(Sex))
data.for.intro.models %>% group_by(StudyID) %>% tally()

## Evolhist
with(data.for.evolhist.models, table(Kingsolver_traits, Sex))
with(data.for.evolhist.models, table(Kingsolver_traits))
with(data.for.evolhist.models, table(Sex))
data.for.evolhist.models %>% group_by(StudyID) %>% tally()

## Determinant Models ----

### 1. Ecology models ----
# Ecology model
## Is there a difference between the slopes? 

### fix structure
data.for.ecology.models$method <- as.factor(data.for.ecology.models$method)

### Remove 'Both' sex category because duplicates
data.for.ecology.models <- data.for.ecology.models %>% filter(Sex %in% c("M", "F"))  

## Model with everything
(ecology.full <- glmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.ecology.models, family = binomial)) %>% summary()

## this is for the means that I report in the text
(one.slope <- data.for.ecology.models %>% filter(method == "south")) %>% summary()
(both.slopes <- data.for.ecology.models %>% filter(method == "all")) %>% summary()


### 2. Intro models ----
## Is there a difference between studies with only natural vs w intro

### Fix strucutre
data.for.intro.models$method <- as.factor(data.for.intro.models$method)

### Remove 'Both' sex category because duplicates
data.for.intro.models <- data.for.intro.models %>% filter(Sex %in% c("M", "F"))  

## Model with everything 
(intro.full <- glmer(R.2 ~ method* + (1|StudyID), data = data.for.intro.models, family = binomial)) %>% summary()

## this is for the means that I report in the text
intro.data.no.colour %>% filter(method == "all") %>% summary()
intro.data.no.colour %>% filter(method == "intro") %>% summary()

## this is for the means that I would report in text (but I don't - "Strict")
only_natural <- data.for.intro.models %>% filter(method == "only_natural")
mean(only_natural$R.2)
sd(only_natural$R.2)
natural_and_intro <- data.for.intro.models %>% filter(method == "all") 
mean(natural_and_intro$R.2)
sd(natural_and_intro$R.2)

### intro BROAD model ----
## Is there a difference between studies with only natural vs w intro

### Fix structure
data.for.intro.models.broad$method <- as.factor(data.for.intro.models.broad$method)

### Remove 'Both' sex category because duplicates
data.for.intro.models.broad <- data.for.intro.models.broad %>% filter(Sex %in% c("M", "F"))  

## Model with everything 
(intro.full.broad <- glmer(R.2 ~ method + (1|StudyID), data = data.for.intro.models.broad, family = binomial)) %>% summary()


## for means/sd reported in text
only_natural_broad <- data.for.intro.models.broad %>% filter(method == "only_natural_broad") 
mean(only_natural_broad$R.2)
sd(only_natural_broad$R.2)
natural_and_intro_broad <- data.for.intro.models.broad %>% filter(method == "all") 
mean(natural_and_intro_broad$R.2)
sd(natural_and_intro_broad$R.2)

### 3. Evolutionary history models ----
## Is there a difference when only w pops in the caroni vs also in the Oropuche? 

## These are all singular fits

### fix structure
data.for.evolhist.models$method <- as.factor(data.for.evolhist.models$method)

### Remove 'Both' sex category because duplicates
data.for.evolhist.models <- data.for.evolhist.models %>% filter(Sex %in% c("M", "F"))  

## Make dataframes to remove colour
evolhist.data.no.colour <- data.for.evolhist.models %>% filter(!Kingsolver_traits == "Colour")

## Model with everything 
(evolhist.full <- glmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.evolhist.models, family = binomial)) %>% summary()

## means/sd reported in text:
one_drainage <- data.for.evolhist.models %>% filter(method == "caroni") 
mean(one_drainage$R.2)
sd(one_drainage$R.2)

both_drainages <- data.for.evolhist.models %>% filter(method == "both.drainages")
mean(both_drainages$R.2)
sd(both_drainages$R.2)


## Troubleshooting ----
### Evolhist are all singular, so here we are comparing our models that are not singular to glms/lmer, to see if it changes anything.

#### (It seems to me like glm/glmer are consistent, so we should be able to use the evolutionary history model and say that we checked (?))

#### Note that I can't plot lmer w an sjplot (cause lmer is plotted as an Estimate, whereas the glm/glmer are odds ratios). I think that the glm/glmer are more consistent anyway (?) but I am only basing this on agreement re significance.

## Compare intro
### This one is not singular (so is an example)
summary(intro.no.colour)
(intro.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = intro.data.no.colour)) %>% summary()
Anova(intro.full.lmer, type = 2)
(intro.full.glm <- glm(R.2 ~ method*Sex, data = intro.data.no.colour, family = binomial)) %>% summary()

## Compare ecology
summary(ecology.no.colour)
(ecology.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = ecology.data.no.colour)) %>% summary()
Anova(ecology.full.lmer, type = 3)
(ecology.full.glm <- glm(R.2 ~ method*Sex, data = ecology.data.no.colour, family = binomial)) %>% summary()


### This is singular (need to decide with this one)

summary(evolhist.full)
(evolhist.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.evolhist.models)) %>% summary()
Anova(evolhist.full.lmer, type = 3)
(evolhist.full.glm <- glm(R.2 ~ method*Sex, data = data.for.evolhist.models, family = binomial)) %>% summary()

car::Anova(evolhist.full.lmer, type = "II")
car::Anova(evolhist.full.glm, type = "II")

summary(evolhist.no.colour)
(evolhist.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = evolhist.data.no.colour)) %>% summary()
Anova(evolhist.full.lmer, type = 3)
(evolhist.full.glm <- glm(R.2 ~ method*Sex, data = evolhist.data.no.colour, family = binomial)) %>% summary()

plot_models(evolhist.no.colour, evolhist.full.glm)


## first figure ----

test <- read.csv("testPG_figure.csv", fileEncoding="UTF-8-BOM")
test$population <- as.factor(test$population)

(highparal <- test %>% filter(facet == 'a') %>% 
    ggplot(aes(x = population, y = mean, color = predation)) +
    geom_point(size = 4) +
    labs(title = "Example of high R2 (0.984)",
         subtitle = #"Study = Reddon et al. 2018, \n
           # "Trait = Male standard length (mm)"
    ) +
    scale_colour_manual(values =  c("black", "grey")) +
    theme_classic() +
    theme_bw() +
    labs(y = "Male standard length (mm)",
         x = "\nPopulation") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 16)))


(lowparal <- test %>% filter(facet == 'b') %>% 
    ggplot(aes(x = population, y = mean, color = predation)) +
    geom_point(size = 4) +
    labs(title = "Example of low R2 (0.000)",
         subtitle = #"Study = Easty et al. 2011, \n
           # "Trait = Mean relative area of black (%)"
    ) +
    scale_colour_manual(values =  c("black", "grey")) +
    theme_classic() +
    theme_bw() +
    labs(y = "Mean relative area of black (%)",
         x = "\nPopulation") +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 16)))

ggarrange(highparal, lowparal, common.legend = TRUE, legend = "bottom")

sup3 <-
  data.all %>% 
  filter(StudyType %in% c("Common Garden (F2)", "Wildcaught")) %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1)) +
  facet_grid(Sex ~ StudyType, scales = "free")


sup2 <-
  data.all %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(angle = 45, vjust= 1, hjust = 1)) +
  facet_grid(Sex ~ Kingsolver_traits, scales = "free")



ggarrange(highparal, lowparal, common.legend = TRUE, legend = "bottom")

highparal <- highparal + annotate("segment", x = 1, xend = 3, y = 12.01, yend = 13.5)
highparal <- highparal + annotate("segment", x = 2, xend = 4, y = 12.11, yend = 13.76)
highparal

lowparal <- lowparal + annotate("segment", x = 1, xend = 3, y = 8.74, yend = 9.86)
lowparal <- lowparal + annotate("segment", x = 2, xend = 4, y = 8.20, yend = 7.0)

lowparal

ggarrange(highparal, lowparal, common.legend = TRUE, legend = "bottom")


# histograms in manuscript ---- 

data.all$overall <- "Overall"
overall_hist <- 
  data.all %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  #ggtitle("Overall (n = 446)") +
  theme_bw() +
  facet_wrap(. ~ overall) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

wildcaught_hist <-
  data.all %>% 
  filter(StudyType == "Wildcaught") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  #ggtitle("Wildcaught (n = 373)") +
  facet_wrap(. ~ StudyType) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

cg_hist <- 
  data.all %>% 
  filter(StudyType == "Common Garden (F2)") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  #ggtitle("Common Garden F1/F2 (n = 73)") +
  facet_wrap(. ~ StudyType) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

m_hist <- 
  data.all %>% 
  filter(Sex == "M") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 #fill="#E6E6E6", 
                 colour = "black", fill = "#E6E6E6", size = 1) +
  labs(#x = expression(paste(R^2)),
    y = "Frequency", 
    # title = "Males (n = 274)",
    #subtitle = "colour traits included"
  ) +
  theme_bw() +
  facet_wrap(. ~ Sex) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

m_hist_no_colour <- 
  data.all.no.colour %>% 
  filter(Sex == "M") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 #fill="#E6E6E6", 
                 colour = "black", fill = "#E6E6E6", size = 1) +
  labs(x = expression(paste(R^2)),
       y = "Frequency",
       #title = "Males (n = 165)",
       #subtitle = "colour traits excluded"
  ) +
  facet_wrap(. ~ Sex) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

f_hist_w_other <-
  data.all %>% 
  filter(Sex == "F") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 #fill="#E6E6E6", 
                 colour = "black", fill = "#E6E6E6", size = 1) +
  #xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  # ggtitle("Females (n = 172)") +
  facet_wrap(. ~ Sex) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

#f_hist_no_other <-
#data.all %>% 
#  filter(Sex == "F" &
#           !Kingsolver_traits == "Other") %>% 
#  ggplot(aes(x = R.2)) +
#  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
#                 #fill="#E6E6E6", 
#                 colour = "black", fill = "#E6E6E6", size = 1) +
#  #xlab(expression(paste(R^2))) +
#  ylab("Frequency") +
#  # ggtitle("Females (n = 172)") +
#  facet_wrap(. ~ Sex) +
#  theme_bw() +
#  theme(
#    strip.text = element_text(size = 13),
#    axis.title = element_text(size = 24))


top.fig.2 <- cowplot::plot_grid("", overall_hist + theme(axis.title = element_blank()),
                                wildcaught_hist + theme(axis.title = element_blank()),
                                cg_hist + theme(axis.title = element_blank()), 
                                nrow = 1, labels = c("", "A", "B", "C"), 
                                rel_widths = c(0.15, 1, 1, 1))

bottom.fig.2 <- cowplot::plot_grid("", m_hist + theme(axis.title = element_blank()),
                                   m_hist_no_colour + theme(axis.title.y = element_blank()),
                                   f_hist_w_other + theme(axis.title = element_blank()),
                                   
                                   #f_hist_no_other + theme(axis.title.y = element_blank()), 
                                   nrow = 1,
                                   labels = c("", "D", "E", "F", "G"),
                                   align = "h", axis = "bt",
                                   rel_widths = c(0.15, 1, 1,1))

bottom.fig.2


fig2 <- cowplot::plot_grid(top.fig.2, bottom.fig.2, nrow = 2, rel_heights = c(1,1))

fig2

fig2 <- fig2 + cowplot::draw_label("Frequency", x=  0, y=0.5, vjust= 1, angle=90, size = 24)

fig2

(fig3<-
    data.all %>% 
    ggplot(aes(x = R.2)) +
    geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                   fill="#E6E6E6", colour = "black", size = 1) +
    xlab(expression(paste(R^2))) +
    ylab("Frequency") +
    theme_bw() +
    facet_wrap(. ~ Kingsolver_traits, nrow = 2, scales = "free") +
    theme(
      strip.text = element_text(size = 13),
      axis.title = element_text(size = 24)))

# table for supplement ----

ok <- data.all %>% group_by(Study)
ok.test <- mutate(ok, min(R.2))
ok.test <- mutate(ok.test, max(R.2))

ok2 <- data.frame(data.all %>% group_by(Study, Poptype_broad) %>% tally())

check <- rbind(ok.test, ok2)

No.traits <- data.all %>% group_by(Study, Trait) %>% tally() %>% 
  transmute(No.traits = sum(n)) %>% 
  distinct(No.traits, .keep_all = FALSE) 
View(No.traits)

Slope.names <- spreadsheet.data %>% group_by(Study, Slope) %>% count()
Slope.names <- data.frame(aggregate(Slope ~ Study, data = Slope.names, paste, collapse = ","))

wat <- left_join(No.traits, Slope.names)

Drainage.names <- spreadsheet.data %>% group_by(Study, Drainage) %>% count()
Drainage.names <- data.frame(aggregate(Drainage ~ Study, data = Drainage.names, paste, collapse = ","))

wat <- left_join(wat, Drainage.names)

Pop.names <- spreadsheet.data %>% group_by(Study, Poptype_broad) %>% count()
Pop.names <- data.frame(aggregate(Poptype_broad ~ Study, data = Pop.names, paste, collapse = ","))

wat <- left_join(wat, Pop.names)


Populations.names <- spreadsheet.data %>% group_by(Study, Poptype_broad, Population) %>% 
  filter(Poptype_broad == "Introduction") %>% count()
Populations.names <- data.frame(aggregate(Population ~ Study, data = Populations.names, paste, collapse = ","))

wat <- left_join(wat, Populations.names)

View(wat)

write.csv(wat, "Supplement_table.csv")

# last fig ----



kk <- list("Male with colour (274)" = Malewithcolour$R.2, "Male without colour (165)" = Malenocolour$R.2, 
           "Female (172)" = Femalesex$R.2, 
           "Colour (109)" = Colourtrait$R.2, "Life history (48)" = LHtrait$R.2, "Size (47)" = Sizetrait$R.2, 
           "Morphology (42)" = Morphtrait$R.2, "Other (3)" = Othertrait$R.2, "Physiology (64)" = Phystrait$R.2,
           "Common Garden (F2) (70)" = cg$R.2, "Wild caught (373)" = wc$R.2, "Within one slope (176)" = one.slope$R.2, "Between both slopes (176)" = both.slopes$R.2, 
           "Only natural (183)" = only_natural$R.2, "Natural and introduced (183)" = natural_and_intro$R.2, 
           "Within one drainage (258)" = one_drainage$R.2, "Between both drainages (258)" = both_drainages$R.2)

mean <- (as.data.frame(sapply(kk, mean)))
mean <- cbind(Factor = rownames(mean), mean)
rownames(mean) <- NULL
colnames(mean)[2] <- "mean"

standev <- (as.data.frame(sapply(kk, sd)))
standev <- cbind(Factor = rownames(standev), standev)
rownames(standev) <- NULL
colnames(standev)[2] <- "sd"

larger_table <- inner_join(mean, standev, by = "Factor")

larger_table$minsd <- larger_table$mean - larger_table$sd
larger_table$maxsd <- larger_table$mean + larger_table$sd

larger_table$order <- c(1,2,3,4,5,6,7,8, 9, 10, 11, 12, 13, 14, 15, 16, 17)
larger_table$Question <- c("Sex", "Sex", "Sex", 
                           "Trait type", "Trait type", "Trait type", "Trait type",
                           "Trait type", "Trait type",
                           "Rearing environment", "Rearing environment",
                           "Ecological complexity", "Ecological complexity",
                           "Time since colonization", "Time since colonization",
                           "Evolutionary history", "Evolutionary history")

larger_table$Question  <- with(larger_table, reorder(Question, rev(order)))

vertical_lines <- c(3.5, 9.5, 11.5, 13.5, 15.5)

(lastfig2b<-
    larger_table %>% ggplot(x = Factor, aes(reorder(Factor, order), y = mean)) +
    geom_point(size = 2) +
    geom_linerange(aes(x = Factor, ymin = minsd, ymax = maxsd), size = 0.5) +
    
    annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf, fill = "yellow", alpha = 0.2) +
    annotate("rect", xmin = 3.5, xmax = 9.5, ymin = -Inf, ymax = Inf, fill = "orange", alpha = .4) +
    annotate("rect", xmin = 9.5, xmax = 11.5, ymin = -Inf, ymax = Inf, fill = "palegreen", alpha = .4) +
    annotate("rect", xmin = 11.5, xmax = 13.5, ymin = -Inf, ymax = Inf, fill = "#0404B7", alpha = .3) +
    annotate("rect", xmin = 13.5, xmax = 15.5, ymin = -Inf, ymax = Inf, fill = "purple", alpha = .4) +
    annotate("rect", xmin = 15.5, xmax = 17.5, ymin = -Inf, ymax = Inf, fill = "hot pink", alpha = .3) +
    
    
    labs(x = "Factor", 
         y = expression(paste(R^2))) +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 24)))

# hist and boxplot figure ----

levels(data.for.intro.models.broad$method)

data.intro.broad <- data.for.intro.models.broad
levels(data.intro.broad$method) <- c("Natural and introduced",
                                     "Only natural")

(intro_hist_broad_a <-
    data.intro.broad %>% 
    filter(method == "Natural and introduced") %>% 
    ggplot(aes(x = R.2)) +
    geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                   fill="#E6E6E6", colour = "black", size = 1) +
    facet_wrap(.~method) +
    labs(x = expression(paste(R^2))) +
    theme_bw() +
    theme(
      axis.ticks.x = element_blank(),
    ) +
    xlim(-0.1,1) +
    theme(axis.title = element_blank(),
          strip.text = element_text(size = 13),
          axis.text = element_text(size = 10)))

intro_hist_broad_a

(intro_hist_broad_b <-
    data.intro.broad %>% 
    filter(method == "Only natural") %>% 
    ggplot(aes(x = R.2)) +
    geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                   fill="#E6E6E6", colour = "black", size = 1) +
    facet_wrap(.~method) +
    labs(x = expression(paste(R^2))) +
    theme_bw() +
    theme(
      axis.ticks.x = element_blank(),
    ) +
    xlim(-0.1,1) +
    #theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
    theme(axis.title = element_blank(),
          strip.text = element_text(size = 13),
          axis.text = element_text(size = 10)))



intro_box <-
  data.for.intro.models.broad %>% 
  ggplot(aes(x = R.2)) +
  geom_boxplot(size = 1) +
  #geom_point(aes(x = method, y = R.2), size = 3, alpha = 0.3) +
  facet_wrap(~method, ncol =1) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  xlim(-0.1,1) +
  theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm"))+
  theme(axis.title = element_blank())
intro_box

intro_box_cow <- cowplot::plot_grid(intro_hist_broad_a, intro_box,
                                    intro_hist_broad_b, nrow = 3, ncol = 1,
                                    rel_heights = c(4,1, 4),
                                    align = "v", axis = "lr", 
                                    labels = c("C", "", "F"))
intro_box_cow

## ecology

levels(data.for.ecology.models$method)

data.eco <- data.for.ecology.models
levels(data.eco$method) <- c("Between both slopes",
                             "Only the southern slope")

(ecology_hist_a <-
    data.eco %>% 
    filter(method == "Between both slopes") %>% 
    ggplot(aes(x = R.2)) +
    geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                   fill="#E6E6E6", colour = "black", size = 1) +
    xlab(expression(paste(R^2))) +
    ylab("Frequency") +
    facet_wrap(~method) +
    labs(x = "") +
    theme_bw() +
    theme(
      axis.ticks.x = element_blank(),
    ) +
    xlim(-0.1,1) +
    # theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
    theme(axis.title = element_blank(),
          strip.text = element_text(size = 13),
          axis.text = element_text(size = 10)))

(ecology_hist_b <-
    data.eco %>% 
    filter(method == "Only the southern slope") %>% 
    ggplot(aes(x = R.2)) +
    geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                   fill="#E6E6E6", colour = "black", size = 1) +
    xlab(expression(paste(R^2))) +
    ylab("Frequency") +
    facet_wrap(~method) +
    labs(x = "") +
    theme_bw() +
    theme(
      axis.ticks.x = element_blank(),
    ) +
    xlim(-0.1,1) +
    # theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
    theme(axis.title = element_blank(),
          strip.text = element_text(size = 13),
          axis.text = element_text(size = 10)))

ecology_box <-
  data.for.ecology.models %>% 
  ggplot(aes(x = R.2)) +
  geom_boxplot(size = 1) +
  #geom_point(aes(x = method, y = R.2), size = 3, alpha = 0.3) +
  facet_wrap(~method, ncol =1) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  xlim(-0.1,1) +
  theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm"))+
  theme(axis.title = element_blank())

ecology_box_cow <- cowplot::plot_grid(ecology_hist_a, ecology_box,
                                      ecology_hist_b, nrow = 3, ncol = 1,
                                      rel_heights = c(4,1, 4),
                                      align = "v", axis = "lr", 
                                      labels = c("A", "", "D"))

ecology_box_cow

## evolutionary history

data.evo <- data.for.evolhist.models
levels(data.evo$method) <- c("Between both drainages",
                             "Only the Caroni drainage")
(evolhist_hist_a <-
    data.evo %>% 
    filter(method == "Between both drainages") %>% 
    ggplot(aes(x = R.2)) +
    geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                   fill="#E6E6E6", colour = "black", size = 1) +
    
    facet_wrap(.~method) +
    labs(x = "") +
    theme_bw() +
    xlab(expression(paste(R^2))) +
    ylab("Frequency") +
    theme(
      axis.title = element_blank(),
      strip.text = element_text(size = 13),
      axis.text = element_text(size = 10)
    ) +
    xlim(-0.1,1) 
  #+  theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
)

(evolhist_hist_b <-
    data.evo %>% 
    filter(method == "Only the Caroni drainage") %>% 
    ggplot(aes(x = R.2)) +
    geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                   fill="#E6E6E6", colour = "black", size = 1) +
    
    facet_wrap(.~method) +
    labs(x = "") +
    theme_bw() +
    xlab(expression(paste(R^2))) +
    ylab("Frequency") +
    theme(
      axis.title = element_blank(),
      strip.text = element_text(size = 13),
      axis.text = element_text(size = 10)
    ) +
    xlim(-0.1,1) 
  #+  theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
)

evolhist_box <-
  data.for.evolhist.models %>% 
  ggplot(aes(x = R.2)) +
  geom_boxplot(size = 1) +
  #geom_point(aes(x = method, y = R.2), size = 3, alpha = 0.3) +
  facet_wrap(~method, ncol =1) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) +
  xlim(-0.1,1) +
  theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm"))+
  theme(axis.title = element_blank())


evolhist_box_cow <- cowplot::plot_grid(evolhist_hist_a, evolhist_box, 
                                       evolhist_hist_b, nrow = 3, ncol = 1,
                                       rel_heights = c(4,1, 4),
                                       align = "v", axis = "lr", 
                                       labels = c("B", "", "E"))

evolhist_box_cow

fig5 <- cowplot::plot_grid(ecology_box_cow, evolhist_box_cow, intro_box_cow, nrow = 1,
                           align = "hv", axis = "tblr",
                           rel_widths = c(1,1,1))
tiff("fig5.pg.tiff", width = 8, height = 4, units = "in", res = 600)
fig5
dev.off()




""







data.all %>% 
  ggplot(aes(x = Sex, y = R.2)) + 
  geom_flat_violin(size = 1, position = position_nudge(x = 0.0, y = 0), adjust = 2) +
  geom_jitter(position = position_nudge(x = - .1, y = 0)) +
  geom_boxplot(size = 1, aes(x = Sex, y = R.2) ,
               alpha = 0.3, width = 0.1) +
  theme_bw() +
  facet_wrap(~Kingsolver_traits) +
  labs(y = expression(paste(R^2)),
       x = "\nMethod") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  theme(legend.position = "none")



data.renamed <- data.all
levels(data.renamed$Kingsolver_traits) <- c("Life history",
                                            "Behaviour",
                                            "Colour",
                                            "Other",
                                            "Morphology",
                                            "Physiology",
                                            "Size")


(figa3 <- 
    data.renamed %>% 
    filter(StudyType %in% c("Wildcaught", "Common Garden (F2)")) %>% 
    filter(Sex %in% c("M", "F")) %>% 
    ggplot(aes(x = R.2)) +
    geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                   fill="#E6E6E6", colour = "black", size = 1) +
    xlab(expression(paste(R^2))) +
    ylab("Frequency") +
    theme_bw() +
    facet_grid(Sex~StudyType, scales = "free") +
    theme(
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, vjust = 0.5),
      axis.title.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      strip.text = element_text(size = 13),
      plot.title = element_text(size = 24)))

tiff("figa3.pg.tiff", res = 600, units = "in", height = 6, width = 8)
figa3
dev.off()

# figure in supplement ----

(figa2 <- 
    data.renamed %>% 
    filter(Sex %in% c("M", "F")) %>% 
    filter(!Kingsolver_traits %in% c("Colour", "Other")) %>% 
    ggplot(aes(x = R.2)) +
    geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                   fill="#E6E6E6", colour = "black", size = 0.75) +
    xlab(expression(paste(R^2))) +
    ylab("Frequency") +
    theme_bw() +
    facet_grid(Sex~Kingsolver_traits, scales = "free") +
    theme(
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, vjust = 0.5),
      axis.title.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      strip.text = element_text(size = 13),
      plot.title = element_text(size = 24)))

tiff("figa2.pg.tiff", res = 600, units = "in", height = 6, width = 8)
figa2
dev.off()

# traits plot ----

behav_hist <-
  data.all %>% 
  filter(Kingsolver_traits == "Behaviour") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  #ggtitle("Overall (n = 446)") +
  theme_bw() +
  facet_wrap(.~ Kingsolver_traits) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_blank())

colour_hist <-
  data.all %>% 
  filter(Kingsolver_traits == "Colour") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  #ggtitle("Overall (n = 446)") +
  theme_bw() +
  facet_wrap(.~ Kingsolver_traits) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_blank())

other_hist <-
  data.all %>% 
  filter(Kingsolver_traits == "Other") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  
  theme_bw() +
  facet_wrap(.~ Kingsolver_traits) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_blank())

data.olh <- data.all
levels(data.olh$Kingsolver_traits) <- c("Life history",
                                        "Behaviour",
                                        "Colour",
                                        "Other",
                                        "Morphology",
                                        "Physiology",
                                        "Size")

levels(data.olh$Kingsolver_traits)

olh_hist <-
  data.olh %>% 
  filter(Kingsolver_traits == "Life history") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  #ggtitle("Overall (n = 446)") +
  theme_bw() +
  facet_wrap(.~ Kingsolver_traits) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_blank())

morpho_hist <-
  data.olh %>% 
  filter(Kingsolver_traits == "Morphology") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  #ggtitle("Overall (n = 446)") +
  theme_bw() +
  facet_wrap(.~ Kingsolver_traits) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_blank())

phys_hist <-
  data.all %>% 
  filter(Kingsolver_traits == "Physiology") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  #ggtitle("Overall (n = 446)") +
  theme_bw() +
  xlab(expression(paste(R^2))) +
  facet_wrap(.~ Kingsolver_traits) +
  theme(
    strip.text = element_text(size = 13),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 24))

size_hist <-
  data.all %>% 
  filter(Kingsolver_traits == "Size") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  #ggtitle("Overall (n = 446)") +
  theme_bw() +
  facet_wrap(.~ Kingsolver_traits) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_blank())


toptraits <- cowplot::plot_grid("", behav_hist, colour_hist, other_hist, olh_hist, nrow = 1, 
                                rel_widths = c(0.15, 1, 1, 1, 1))
bottomtraits <- cowplot::plot_grid("", "", morpho_hist, phys_hist, size_hist, "", nrow = 1,
                                   rel_widths = c(0.15, 0.4, 0.8, 0.8, 0.8, 0.4), 
                                   align = "h", axis = "bt")
traitsplot <- cowplot::plot_grid(toptraits,bottomtraits, nrow = 2, align = "bt", rel_heights = c(0.9, 1))
traitsplot <- traitsplot + cowplot::draw_label("Frequency", x=  0, y=0.5, vjust= 1, angle=90, size = 24)

tiff("fig4.pg.tiff", res= 600, height = 6, width = 8, units = "in")
traitsplot
dev.off()

##

##



## OLD AGAIN ----

## updates 2021-10-21 as I'm cleaning AH

## trait type plots ----
## here is a general plot with colour
data.all.traits %>% 
  ggplot(aes(x = Kingsolver_traits, y = R.2, fill = Kingsolver_traits)) + 
  geom_jitter(aes(color = Kingsolver_traits), width = .1) +
  geom_boxplot(size = 1, aes(x = Kingsolver_traits, y = R.2) ,
               alpha = 0.3, width = 0.3) +
  scale_colour_manual(values =  c("#FFB93C", "#457111","#15899A", "#BD8DC3", "#D07D7D", "dark gray", "pink")) +
  scale_fill_manual(values =  c("#FFB93C", "#457111","#15899A", "#BD8DC3", "#D07D7D", "dark gray", "pink")) +
  theme_classic() 

data.all.rear %>% 
  ggplot(aes(x = StudyType, y = R.2, color = Kingsolver_traits)) +
  geom_boxplot() + 
  geom_jitter()

data.all.no.colour %>% 
  ggplot(aes(x = Kingsolver_traits, y = R.2, fill = Kingsolver_traits)) + 
  geom_jitter(aes(color = Kingsolver_traits), width = .1) +
  geom_boxplot(size = 1, aes(x = Kingsolver_traits, y = R.2) ,
               alpha = 0.3, width = 0.3) +
  scale_colour_manual(values =  c("#FFB93C", "#457111","#15899A", "#BD8DC3", "#D07D7D", "dark gray")) +
  scale_fill_manual(values =  c("#FFB93C", "#457111","#15899A", "#BD8DC3", "#D07D7D", "dark gray")) +
  theme_classic() 

## sex plots ----
data.all.traits %>% 
  ggplot(aes(x = Sex, y = R.2, fill = Sex)) + 
  geom_jitter(aes(color = Sex), width = 0.1, size = 2) +
  geom_boxplot(size = 1, aes(x = Sex, y = R.2), alpha = 0.3) +
  scale_colour_manual(values =  c("#FFB93C", "#457111")) +
  scale_fill_manual(values =  c("#FFB93C", "#457111")) +
  theme_classic()  

data.all.no.colour %>% 
  ggplot(aes(x = Sex, y = R.2, fill = Sex)) + 
  geom_jitter(aes(color = Sex), width = 0.1, size = 2) +
  geom_boxplot(size = 1, aes(x = Sex, y = R.2), alpha = 0.3) +
  scale_colour_manual(values =  c("#FFB93C", "#457111")) +
  scale_fill_manual(values =  c("#FFB93C", "#457111")) +
  theme_classic()  

data.all %>% filter(Kingsolver_traits %in% c("Behaviour", "Colour", "Other_life_history", "Size") &
                      StudyType %in% c("Common Garden (F2)", "Wildcaught"))  %>% 
  ggplot(aes(x = StudyType, y = R.2)) +
  geom_boxplot() +
  geom_jitter(width = .1, size = 3, alpha = 0.5, aes(colour = Kingsolver_traits)) +
  facet_wrap(~Kingsolver_traits) +
  theme_bw()


## quesiton-specific plots
#### Plots ----

(ecology.full.plot <-
   ecology.data.no.colour %>% 
   ggplot(aes(x = method, y = R.2, fill = method)) + 
   geom_flat_violin(size = 1, position = position_nudge(x = 0.0, y = 0), adjust = 2) +
   geom_jitter(aes(color = method), position = position_nudge(x = - .1, y = 0)) +
   geom_boxplot(size = 1, aes(x = method, y = R.2) ,
                alpha = 0.3, width = 0.1) +
   scale_colour_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
   scale_fill_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
   theme_bw()+
   theme(legend.position = "none"))


(ecology.sex.plot <-  
    ecology.data.no.colour %>% 
    ggplot(aes(x = method, y = R.2, fill = method)) + 
    geom_flat_violin(size = 1, position = position_nudge(x = 0.0, y = 0), adjust = 2) +
    geom_jitter(aes(color = method), position = position_nudge(x = - .1, y = 0)) +
    geom_boxplot(size = 1, aes(x = method, y = R.2) ,
                 alpha = 0.3, width = 0.1) +
    scale_colour_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    scale_fill_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    facet_wrap(~Sex) +
    theme_bw() )

ggarrange(ecology.full.plot, ecology.sex.plot, ncol = 1, common.legend= TRUE)



#### Intro plots ----
(intro.full.plot <-
   intro.data.no.colour %>% 
   ggplot(aes(x = method, y = R.2, fill = method)) + 
   geom_flat_violin(size = 1, position = position_nudge(x = 0.0, y = 0), adjust = 2) +
   geom_jitter(aes(color = method), position = position_nudge(x = - .1, y = 0)) +
   geom_boxplot(size = 1, aes(x = method, y = R.2) ,
                alpha = 0.3, width = 0.1) +
   scale_colour_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
   scale_fill_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
   theme_bw()+
   theme(legend.position = "none"))


(intro.sex.plot <-  
    intro.data.no.colour %>% 
    ggplot(aes(x = method, y = R.2, fill = method)) + 
    geom_flat_violin(size = 1, position = position_nudge(x = 0.0, y = 0), adjust = 2) +
    geom_jitter(aes(color = method), position = position_nudge(x = - .1, y = 0)) +
    geom_boxplot(size = 1, aes(x = method, y = R.2) ,
                 alpha = 0.3, width = 0.1) +
    scale_colour_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    scale_fill_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    facet_wrap(~Sex) +
    theme_bw())

ggarrange(intro.full.plot, intro.sex.plot, ncol = 1, common.legend= TRUE)
#

#### intro BROAD plots ----
plot_models(intro.full, intro.no.colour, vline.color = "grey")

intro.data.no.colour.broad %>% 
  ggplot(aes(x = method, y = R.2, fill = Sex)) +
  geom_boxplot(size = 1.5) +
  geom_jitter(size = 3, alpha = 0.1, width = 0.2) +
  #scale_fill_manual(values = c("gray48", "gray90")) +
  facet_wrap(~Sex) +
  theme_bw()

(intro.full.plot <-
    intro.data.no.colour.broad %>% 
    ggplot(aes(x = method, y = R.2, fill = method)) + 
    geom_flat_violin(size = 1, position = position_nudge(x = 0.0, y = 0), adjust = 2) +
    geom_jitter(aes(color = method), position = position_nudge(x = - .1, y = 0)) +
    geom_boxplot(size = 1, aes(x = method, y = R.2) ,
                 alpha = 0.3, width = 0.1) +
    scale_colour_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    scale_fill_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    theme_bw()+
    theme(legend.position = "none"))


(intro.sex.plot <-  
    intro.data.no.colour.broad %>% 
    ggplot(aes(x = method, y = R.2, fill = method)) + 
    geom_flat_violin(size = 1, position = position_nudge(x = 0.0, y = 0), adjust = 2) +
    geom_jitter(aes(color = method), position = position_nudge(x = - .1, y = 0)) +
    geom_boxplot(size = 1, aes(x = method, y = R.2) ,
                 alpha = 0.3, width = 0.1) +
    scale_colour_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    scale_fill_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    facet_wrap(~Sex) +
    theme_bw())

ggarrange(intro.full.plot, intro.sex.plot, ncol = 1, common.legend= TRUE)



#### evolhist plots ----
plot_models(evolhist.full, evolhist.no.colour, vline.color = "grey")


(evolhist.full.plot <-
    evolhist.data.no.colour %>% 
    ggplot(aes(x = method, y = R.2, fill = method)) + 
    geom_flat_violin(size = 1, position = position_nudge(x = 0.0, y = 0), adjust = 2) +
    geom_jitter(aes(color = method), position = position_nudge(x = - .1, y = 0)) +
    geom_boxplot(size = 1, aes(x = method, y = R.2) ,
                 alpha = 0.3, width = 0.1) +
    scale_colour_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    scale_fill_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    theme_bw()+
    theme(legend.position = "none"))


(evolhist.sex.plot <-  
    evolhist.data.no.colour %>% 
    ggplot(aes(x = method, y = R.2, fill = method)) + 
    geom_flat_violin(size = 1, position = position_nudge(x = 0.0, y = 0), adjust = 2) +
    geom_jitter(aes(color = method), position = position_nudge(x = - .1, y = 0)) +
    geom_boxplot(size = 1, aes(x = method, y = R.2) ,
                 alpha = 0.3, width = 0.1) +
    scale_colour_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    scale_fill_manual(values =  c("#15899A", "#BD8DC3", "#E27474")) +
    facet_wrap(~Sex) +
    theme_bw() )

ggarrange(evolhist.full.plot, evolhist.sex.plot, ncol = 1, common.legend= TRUE)



## older than 2021-10-21
#### merge spreadsheet data w BOTH of these R2 to get an ecology data spreadsheet

## filter by sex (because duplicates in 'Both')
data.for.models <- filter(data.for.models, Sex %in% c("M", "F"))

# rename across slope set - because data.for.models will be used below too
data.for.models.across <- data.for.models

## QUESTION 2 EVOUTIONARY HISTORY

## first, make a data frame for caroni only
R2.data.caroni$TraitID <- as.factor(R2.data.caroni$TraitID)
R2.data.among.drainage$TraitID <- as.factor(R2.data.among.drainage$TraitID)

data.for.models.caroni <- left_join(spreadsheet.data, R2.data.caroni,  by = "TraitID")

## filter by sex to remove duplicates from ('Both')
data.for.models.caroni <- filter(data.for.models.caroni, Sex %in% c("M", "F"))

## caroni only (again, idk if this really matters cuz ANOVA was run w only caroni and TraitID not specific 2 drainage)
data.for.models.caroni <- data.for.models.caroni %>% filter(Drainage == "Caroni")

## remove duplicates so only 1 R2 per traitID
data.for.models.caroni <- data.for.models.caroni[!duplicated(data.for.models.caroni$TraitID),]

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
