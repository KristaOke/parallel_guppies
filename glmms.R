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
R2.data.intro.broad <- read.csv(paste(wd, '/Data/TraitR2_intro_broad.csv', sep = ""), header = TRUE, sep = "")


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

## combine R2 and spreadsheet data
#prep for when we bind them together, this variable will go in the model to indicate the type of R2
R2.data.among$method <- "all"
R2.data.south$method <- "south"
R2.data.intro$method <- "only_natural"
R2.data.intro.broad$method <- "only_natural_broad"
R2.data.caroni$method <- "caroni"
R2.data.among.drainage$method <- "both.drainages"

#get relevant data with one entry for each Trait
##I was inclusive with columns, many probably aren't needed so you can cut them out
data.all <- inner_join(spreadsheet.data, R2.data.among, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:52)%>% #I selected only the columns that have information that applies at the trait level
  distinct(TraitID, .keep_all = TRUE) #removes replicated TraitIDs and retains the columns

#repeat for south only R2 and others...
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

#then we can bind them together as we wish! First ecology...
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

#evolutionary history question
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

# Intro, I'm pretty sure these are the right csvs for this question
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


# Intro, I'm pretty sure these are the right csvs for this question
data.for.intro.models.broad<-rbind(data.all,data.intro.broad) %>% 
  arrange(TraitID) %>% #puts them in a nice order
  group_by(TraitID) %>% #groups them for the count
  filter(n() > 1) #filters only trait IDs that have more than 1 entry within the group (n = 204 traits)

# as a note and fail safe I would run the grouping with TraitID on the sex specific dfs to make sure theres two traitID entries for each

## Overall models

## remove 'both'
### (Because when calculated by hand, I would do M/F/Both with same data)
data.all <- data.all %>% filter(Sex %in% c("M", "F"))  

sd(data.all$R.2)
#sd(data.all.no.colour$R.2)

## remove other (because model below will not converge with other)
data.all.traits <- data.all %>% filter(!Kingsolver_traits == "Other")

## creating a dataset w no colour (to compare w those w colour)
data.all.no.colour <- data.all %>% filter(!Kingsolver_traits == "Colour")

## basic histograms
data.all.traits %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1.5) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  theme_bw()

## basic histograms - no colour
data.all.no.colour %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1.5) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  theme_bw()

# colour traits included
w.colour.M <- data.all %>% filter(Sex == "M")
mean(w.colour.M$R.2)
sd(w.colour.M$R.2)

w.colour.F <- data.all %>% filter(Sex == "F")
mean(w.colour.F$R.2)
sd(w.colour.F$R.2)

# colour traits excluded
no.colour.M <- data.all.no.colour %>% filter(Sex == "M")
mean(no.colour.M$R.2)
sd(no.colour.M$R.2)

no.colour.F <- data.all.no.colour %>% filter(Sex == "F")
mean(no.colour.F$R.2)
sd(no.colour.F$R.2)

## overall traits model

## this WILL NOT RUN with other included
(all.model.traits <- glmer(R.2 ~ Kingsolver_traits +  (1|StudyID), data = data.all.traits, family = binomial)) %>% summary()

data.all.traits.rear <- data.all.traits %>% filter(StudyType %in% c("Common Garden (F2)", "Wildcaught"))
(all.model.rearing <- glmer(R.2 ~ StudyType +  (1|StudyID), data = data.all.traits.rear, family = binomial)) %>% summary()

r.squaredGLMM(all.model.traits)

## overall sex model
(all.model.sex <- glmer(R.2 ~ Sex + (1|StudyID), data = data.all, family = binomial)) %>% summary()

r.squaredGLMM(all.model.sex)

confint(all.model.sex)

### Removed colour 
(all.model.sex.no.colour <- glmer(R.2 ~ Sex + (1|StudyID), data = data.all.no.colour, family = binomial)) %>% summary()

r.squaredGLMM(all.model.sex.no.colour)

confint(all.model.sex.no.colour)

## trait type plots
## here is a general plot with colour
data.all.traits %>% 
  ggplot(aes(x = Kingsolver_traits, y = R.2, fill = Kingsolver_traits)) + 
  geom_jitter(aes(color = Kingsolver_traits), width = .1) +
  geom_boxplot(size = 1, aes(x = Kingsolver_traits, y = R.2) ,
               alpha = 0.3, width = 0.3) +
  scale_colour_manual(values =  c("#FFB93C", "#457111","#15899A", "#BD8DC3", "#D07D7D", "dark gray", "pink")) +
  scale_fill_manual(values =  c("#FFB93C", "#457111","#15899A", "#BD8DC3", "#D07D7D", "dark gray", "pink")) +
  theme_classic() 

data.all.traits.rear %>% 
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

## sex plots
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

# Traits

Colourtrait <- data.all %>% filter(Kingsolver_traits == "Colour")
mean(Colourtrait$R.2)
sd(Colourtrait$R.2)

LHtrait <- data.all %>% filter(Kingsolver_traits == "Other_life_history")
mean(LHtrait$R.2)
sd(LHtrait$R.2)

Sizetrait <- data.all %>% filter(Kingsolver_traits == "Size")
mean(Sizetrait$R.2)
sd(Sizetrait$R.2)

Morphtrait <- data.all %>% filter(Kingsolver_traits == "Other_morphology")
mean(Morphtrait$R.2)
sd(Morphtrait$R.2)

Othertrait <- data.all %>% filter(Kingsolver_traits == "Other")
mean(Othertrait$R.2)
sd(Othertrait$R.2)

# Rearing enviro
w.colour.cg <- data.all %>% filter(StudyType == "Common Garden (F2)")
mean(w.colour.cg$R.2)
sd(w.colour.cg$R.2)

w.colour.wc <- data.all %>% filter(StudyType == "Wildcaught")
mean(w.colour.wc$R.2)
sd(w.colour.wc$R.2)

no.colour.cg <- data.all.no.colour %>% filter(StudyType == "Common Garden (F2)")
mean(no.colour.cg$R.2)
sd(no.colour.cg$R.2)

no.colour.wc <- data.all.no.colour %>% filter(StudyType == "Wildcaught")
mean(no.colour.wc$R.2)
sd(no.colour.wc$R.2)


## Determinants

### Sample size tables

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

### Models

### 1. Ecology models 
# Ecology model
## Is there a difference between the slopes? 

### Remove 'Both' sex category because duplicates
data.for.ecology.models <- data.for.ecology.models %>% filter(Sex %in% c("M", "F"))  

## Remove other (to be consistent with above, but idk)
#data.for.ecology.models %>% filter(!Kingsolver_traits == "Other")

## Make dataframes to remove colour
ecology.data.no.colour <- data.for.ecology.models %>% filter(!Kingsolver_traits == "Colour")

## Model with everything
(ecology.full <- glmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.ecology.models, family = binomial)) %>% summary()

r.squaredGLMM(ecology.full)

### Model without colour
(ecology.no.colour <- glmer(R.2 ~ method*Sex + (1|StudyID), data = ecology.data.no.colour, family = binomial)) %>% summary()

r.squaredGLMM(ecology.no.colour)


## effect size plot with all models
plot_models(ecology.full, ecology.no.colour, vline.color = "grey")

## this is for the means that I report in the text
# colour traits included
w.colour.south <- data.for.ecology.models %>% filter(method == "south")
mean(w.colour.south$R.2)
sd(w.colour.south$R.2)
w.colour.all <- data.for.ecology.models %>% filter(method == "all") 
mean(w.colour.all$R.2)
sd(w.colour.all$R.2)

# colour traits excluded
no.colour.south <- ecology.data.no.colour %>% filter(method == "south")
mean(no.colour.south$R.2)
sd(no.colour.south$R.2)
no.colour.all <- ecology.data.no.colour %>% filter(method == "all") 
mean(no.colour.all$R.2)
sd(no.colour.all$R.2)


## Plots

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



### 2. Intro models
# intro model
## Is there a difference between studies with only natural vs w intro

### Remove 'Both' sex category because duplicates
data.for.intro.models <- data.for.intro.models %>% filter(Sex %in% c("M", "F"))  

### Remove 'other' to be consistent
#data.for.intro.models <- data.for.intro.models %>% filter(!Kingsolver_traits == "Other")

## Make dataframes to remove colour
intro.data.no.colour <- data.for.intro.models %>% filter(!Kingsolver_traits == "Colour")

## Model with everything 
(intro.full <- glmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.intro.models, family = binomial)) %>% summary()

r.squaredGLMM(intro.full)

### Model without colour
(intro.no.colour <- glmer(R.2 ~ method*Sex + (1|StudyID), data = intro.data.no.colour, family = binomial)) %>% summary()

r.squaredGLMM(intro.no.colour)


## effect size plot with all models
plot_models(intro.full, intro.no.colour, vline.color = "grey")

intro.data.no.colour %>% 
  ggplot(aes(x = method, y = R.2)) +
  geom_boxplot(size = 1.5) +
  geom_jitter(size = 3, alpha = 0.1, width = 0.2) +
  scale_fill_manual(values = c("gray48", "gray90")) +
  #facet_wrap(~Sex) +
  theme_bw()


## this is for the means that I report in the text
intro.data.no.colour %>% filter(method == "all") %>% summary()
intro.data.no.colour %>% filter(method == "intro") %>% summary()

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

## this is for the means that I report in the text
# colour traits included
w.colour.nat <- data.for.intro.models %>% filter(method == "only_natural")
mean(w.colour.nat$R.2)
sd(w.colour.nat$R.2)
w.colour.all <- data.for.intro.models %>% filter(method == "all") 
mean(w.colour.all$R.2)
sd(w.colour.all$R.2)

# colour traits excluded
no.colour.nat <- intro.data.no.colour %>% filter(method == "only_natural")
mean(no.colour.nat$R.2)
sd(no.colour.nat$R.2)
no.colour.all <- intro.data.no.colour %>% filter(method == "all") 
mean(no.colour.all$R.2)
sd(no.colour.all$R.2)


# intro %>% model
## Is there a difference between studies with only natural vs w intro

### Remove 'Both' sex category because duplicates
data.for.intro.models.broad <- data.for.intro.models.broad %>% filter(Sex %in% c("M", "F"))  

### Remove 'other' to be consistent
#data.for.intro.models.broad <- data.for.intro.models.broad %>% filter(!Kingsolver_traits == "Other")

## Make dataframes to remove colour
intro.data.no.colour.broad <- data.for.intro.models.broad %>% filter(!Kingsolver_traits == "Colour")

## Model with everything 
(intro.full.broad <- glmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.intro.models.broad, family = binomial)) %>% summary()

r.squaredGLMM(intro.full.broad)

### Model without colour
(intro.no.colour <- glmer(R.2 ~ method*Sex + (1|StudyID), data = intro.data.no.colour.broad, family = binomial)) %>% summary()

r.squaredGLMM(intro.no.colour)


## effect size plot with all models
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


w.colour.nat.broad <- data.for.intro.models.broad %>% filter(method == "only_natural_broad") 
mean(w.colour.nat.broad$R.2)
sd(w.colour.nat.broad$R.2)
w.colour.all.broad <- data.for.intro.models.broad %>% filter(method == "all") 
mean(w.colour.all.broad$R.2)
sd(w.colour.all.broad$R.2)

# colour traits excluded
no.colour.nat.broad <- intro.data.no.colour.broad %>% filter(method == "only_natural_broad")
mean(no.colour.nat.broad$R.2, na.rm = TRUE)
sd(no.colour.nat.broad$R.2, na.rm = TRUE)
no.colour.all.broad <- intro.data.no.colour.broad %>% filter(method == "all") 
mean(no.colour.all.broad$R.2)
sd(no.colour.all.broad$R.2)


# 3. Evolutionary history models
# evolhist model
## Is there a difference when only w pops in the caroni vs also in the Oropuche? 

## These are all singular fits

### Remove 'Both' sex category because duplicates
data.for.evolhist.models <- data.for.evolhist.models %>% filter(Sex %in% c("M", "F"))  

### Being consistent

#data.for.evolhist.models <- data.for.evolhist.models %>% filter(!Kingsolver_traits == "Other")

## Make dataframes to remove colour
evolhist.data.no.colour <- data.for.evolhist.models %>% filter(!Kingsolver_traits == "Colour")

## Model with everything 
(evolhist.full <- glmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.evolhist.models, family = binomial)) %>% summary()

r.squaredGLMM(evolhist.full)

### Model without colour
(evolhist.no.colour <- glmer(R.2 ~ method*Sex + (1|StudyID), data = evolhist.data.no.colour, family = binomial)) %>% summary()

r.squaredGLMM(evolhist.no.colour)


## effect size plot with all models
plot_models(evolhist.full, evolhist.no.colour, vline.color = "grey")

evolhist.data.no.colour %>% filter(method == "both.drainages") %>% summary()
evolhist.data.no.colour %>% filter(method == "caroni") %>% summary()


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


ok <- evolhist.data.no.colour %>% filter(method == "caroni") 
mean(ok$R.2)
sd(ok$R.2)
ok2 <- evolhist.data.no.colour %>% filter(method == "both.drainages")
mean(ok2$R.2)
sd(ok2$R.2)

ok <- data.for.evolhist.models %>% filter(method == "caroni") 
mean(ok$R.2)
sd(ok$R.2)
ok2 <- data.for.evolhist.models %>% filter(method == "both.drainages")
mean(ok2$R.2)
sd(ok2$R.2)

  


## Troubleshooting 
### Evolhist are all singular, so here we are comparing our models that are not singular to glms/lmer, to see if it changes anything.

#### (It seems to me like glm/glmer are consistent, so we should be able to use the evolutionary history model and say that we checked (?))

#### Note that I can't plot lmer w an sjplot (cause lmer is plotted as an Estimate, whereas the glm/glmer are odds ratios). I think that the glm/glmer are more consistent anyway (?) but I am only basing this on agreement re significance.

## Compare intro
### This one is not singular (so is an example)
summary(intro.no.colour)
(intro.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = intro.data.no.colour)) %>% summary()
Anova(intro.full.lmer, type = 2)
(intro.full.glm <- glm(R.2 ~ method*Sex, data = intro.data.no.colour, family = binomial)) %>% summary()

## Plot models
plot_models(intro.no.colour, intro.full.glm)

## Compare ecology
summary(ecology.no.colour)
(ecology.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = ecology.data.no.colour)) %>% summary()
Anova(ecology.full.lmer, type = 3)
(ecology.full.glm <- glm(R.2 ~ method*Sex, data = ecology.data.no.colour, family = binomial)) %>% summary()

## Plot model to compare effect sizes
plot_models(ecology.no.colour, ecology.full.glm)

#### GLMs seem to be consistent (in terms of p-values and effect sizes) with our glmer that are not singular fits

### This is singular (need to decide with this one)

summary(evolhist.full)
(evolhist.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.evolhist.models)) %>% summary()
Anova(evolhist.full.lmer, type = 3)
(evolhist.full.glm <- glm(R.2 ~ method*Sex, data = data.for.evolhist.models, family = binomial)) %>% summary()

plot_models(evolhist.no.colour, evolhist.full.glm)

summary(evolhist.no.colour)
(evolhist.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = evolhist.data.no.colour)) %>% summary()
Anova(evolhist.full.lmer, type = 3)
(evolhist.full.glm <- glm(R.2 ~ method*Sex, data = evolhist.data.no.colour, family = binomial)) %>% summary()

plot_models(evolhist.no.colour, evolhist.full.glm)


## Plot models (is there a better way to compare effect sizes?)
### (Note that I'm not sure if lmer looks like this because it's bad, or because it would typically be plotted as an estimate, not an odds ratio)
plot_models(ecology.full, ecology.full.glm)

test <- read.csv("testPG_figure.csv", fileEncoding="UTF-8-BOM")
test$population <- as.factor(test$population)

highparal <- test %>% filter(facet == 'a') %>% 
  ggplot(aes(x = population, y = mean, color = predation)) +
  geom_point(size = 4) +
  labs(title = "Example of high R2 (0.984)",
       subtitle = "Study = Reddon et al. 2018, \nTrait = Standard length (mm) (males only)") +
  scale_colour_manual(values =  c("#9F0620", "#69A9C7")) +
  theme_classic()

lowparal <- test %>% filter(facet == 'b') %>% 
  ggplot(aes(x = population, y = mean, color = predation)) +
  geom_point(size = 4) +
  labs(title = "Example of low R2 (0.000)",
       subtitle = "Study = Eastya et al. 2011, \nTrait = Mean relative area of black") +
  scale_colour_manual(values =  c("#9F0620", "#69A9C7")) +
  theme_classic()

ggarrange(highparal, lowparal, common.legend = TRUE, legend = "bottom")


##

##

##



## OLD AGAIN ----
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
