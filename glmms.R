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

## Note: 2021-05-04: Study 41: Female response to male from ancestral population is duplicated. Delete these four rows.

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
library(ggpubr)
library(PupillometryR)

# 2021-06-10 - I have updated this to be consistent with everything on the doc!!

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
R2.data.among$TraitID <- as.factor(R2.data.among$TraitID)
R2.data.south$TraitID <- as.factor(R2.data.south$TraitID)
R2.data.intro$TraitID <- as.factor(R2.data.intro$TraitID)
R2.data.caroni$TraitID <- as.factor(R2.data.caroni$TraitID)
R2.data.among.drainage$TraitID <- as.factor(R2.data.among.drainage$TraitID)

# adding extra columns (broader random)

## these are all of the authors
spreadsheet.data <- 
  spreadsheet.data %>% 
  mutate(first_author = case_when(
    Study == "O'Steen et al. 2002" ~ "O'Steen",  
    Study == "Weese et al. 2010" ~ "Weese",  
    Study == "Burns et al. 2009" ~ "Rodd",  ## this one breask the rules because Burns appears twice
    Study == "Reznick et al. 2005" ~ "Reznick",
    Study == "Evans et al. 2011" ~ "Evans",
    Study == "Kam et al. 2014" ~ "Kam",  # This is a typo in spreadsheet, actually Sandkam
    Study == "Burns and Rodd 2008" ~ "Rodd",  ## burns appears twice
    Study == "Harris et al. 2010" ~ "Harris",
    Study == "Neff et al. 2008" ~ "Neff", 
    Study == "Herbert-Read et al. 2019" ~ "Herbert-Read",
    Study == "Stephenson et al. 2015" ~ "Stephenson",
    Study == "Edenbrow et al. 2017" ~ "Croft",
    Study == "Evans et al. 2003" ~ "Evans",
    Study == "Piyapong et al. 2011" ~ "Piyapong",
    Study == "Magurran and Seghers 1994" ~ "Magurran",
    Study == "Martin and Johnsen 2017" ~ "Martin",
    Study == "Gotanda et al. 2013" ~ "Gotanda",  ## 
    Study == "Hain et al. 2016" ~ "Neff",
    Study == "de Lira et al. inprep" ~ "Hendry",
    Study == "Eastya et al. 2011" ~ "Hendry",
    Study == "Gordon et al. 2017b" ~ "Gordon",
    Study == "Fraser et al. 2015" ~ "Fraser",
    Study == "Devigili et al. 2019" ~ "Devigili",
    Study == "Perez-Jvostov et al. 2012" ~ "Perez-Jvostov",
    Study == "Dial et al. 2016" ~ "Dial",
    Study == "Egset et al. 2011" ~ "Egset",
    Study == "El-Sabaawi et al. 2012" ~ "El-Sabaawi",
    Study == "Huizinga et al. 2009" ~ "Reznick",
    Study == "Reznick et al. 2004" ~ "Reznick",
    Study == "Herbert-Read et al. 2017" ~ "Herbert-Read",
    Study == "Schwartz and Hendry 2007" ~ "Hendry",
    Study == "Zandona et al. 2017" ~ "Zandona",
    Study == "Millar and Hendry 2011" ~ "Hendry",
    Study == "Ioannou et al. 2017" ~ "Ioannou",
    Study == "Reznick and Endler 1982" ~ "Reznick",
    Study == "Reddon et al. 2018" ~ "Reddon",
    Study == "Fischer et al. 2013" ~ "Fischer",
    Study == "Elgee et al. 2010" ~ "Elgee",
    Study == "Gordon et al. 2012" ~ "Gordon",
    Study == "Evans and Magurran 1999" ~ "Evans",
    Study == "Dial et al. 2017" ~ "Dial",
    Study == "Bohr Brask et al. 2019" ~ "Brask",
    Study == "Fischer et al. 2016" ~ "Fischer",
    Study == "Valvo et al. 2019" ~ "Rodd",
    Study == "Reznick and Bryant 2007" ~ "Reznick",
    Study == "Magurran and Seghers" ~ "Magurran",
    Study == "Auer et al. 2018" ~ "Reznick",
    Study == "Croft et al. 2009" ~ "Croft",
    Study == "Zandona et al. 2015" ~ "Zandona"))

## Each of these authors appears more than once (the ones that are # are only represented once)
spreadsheet.data <- 
  spreadsheet.data %>% 
  mutate(first_author_NS = case_when(
    #Study == "O'Steen et al. 2002" ~ "Other",  
    #Study == "Weese et al. 2010" ~ "Other",  
    Study == "Burns et al. 2009" ~ "Rodd",  ## this one breask the rules because Burns appears twice
    Study == "Reznick et al. 2005" ~ "Reznick",
    Study == "Evans et al. 2011" ~ "Evans",
    #Study == "Kam et al. 2014" ~ "Other",  # This is a typo in spreadsheet, actually Sandkam
    Study == "Burns and Rodd 2008" ~ "Rodd",  ## burns appears twice
    #Study == "Harris et al. 2010" ~ "Other",
    Study == "Neff et al. 2008" ~ "Neff", 
    Study == "Herbert-Read et al. 2019" ~ "Herbert-Read",
    #Study == "Stephenson et al. 2015" ~ "Other",
    Study == "Edenbrow et al. 2017" ~ "Croft",
    Study == "Evans et al. 2003" ~ "Evans",
    #Study == "Piyapong et al. 2011" ~ "Other",
    Study == "Magurran and Seghers 1994" ~ "Magurran",
    #Study == "Martin and Johnsen 2017" ~ "Other",
    Study == "Gotanda et al. 2013" ~ "Gotanda",  ## 
    Study == "Hain et al. 2016" ~ "Neff",
    Study == "de Lira et al. inprep" ~ "Hendry",
    Study == "Eastya et al. 2011" ~ "Hendry",
    Study == "Gordon et al. 2017b" ~ "Gordon",
    #Study == "Fraser et al. 2015" ~ "Other",
    #Study == "Devigili et al. 2019" ~ "Other",
    #Study == "Perez-Jvostov et al. 2012" ~ "Other",
    Study == "Dial et al. 2016" ~ "Dial",
    #Study == "Egset et al. 2011" ~ "Other",
    #Study == "El-Sabaawi et al. 2012" ~ "Other",
    Study == "Huizinga et al. 2009" ~ "Reznick",
    Study == "Reznick et al. 2004" ~ "Reznick",
    Study == "Herbert-Read et al. 2017" ~ "Herbert-Read",
    Study == "Schwartz and Hendry 2007" ~ "Hendry",
    Study == "Zandona et al. 2017" ~ "Zandona",
    Study == "Millar and Hendry 2011" ~ "Hendry",
    #Study == "Ioannou et al. 2017" ~ "Other",
    Study == "Reznick and Endler 1982" ~ "Reznick",
    #Study == "Reddon et al. 2018" ~ "Other",
    Study == "Fischer et al. 2013" ~ "Fischer",
    #Study == "Elgee et al. 2010" ~ "Other",
    Study == "Gordon et al. 2012" ~ "Gordon",
    Study == "Evans and Magurran 1999" ~ "Evans",
    Study == "Dial et al. 2017" ~ "Dial",
    #Study == "Bohr Brask et al. 2019" ~ "Other",
    Study == "Fischer et al. 2016" ~ "Fischer",
    Study == "Valvo et al. 2019" ~ "Rodd",
    Study == "Reznick and Bryant 2007" ~ "Reznick",
    Study == "Magurran and Seghers" ~ "Magurran",
    Study == "Auer et al. 2018" ~ "Reznick",
    Study == "Croft et al. 2009" ~ "Croft",
    Study == "Zandona et al. 2015" ~ "Zandona"))

spreadsheet.data <- 
  spreadsheet.data %>% 
  mutate(institution = case_when(
    Study == "O'Steen et al. 2002" ~ "UC Irvine",  
    Study == "Weese et al. 2010" ~ "Maine",  
    Study == "Burns et al. 2009" ~ "Toronto",  
    Study == "Reznick et al. 2005" ~ "UC Riverside",
    Study == "Evans et al. 2011" ~ "Western Australia",
    Study == "Kam et al. 2014" ~ "Simon Fraser",  # This is a typo in spreadsheet, actually Sandkam
    Study == "Burns and Rodd 2008" ~ "Toronto",  
    Study == "Harris et al. 2010" ~ "Lund",
    Study == "Neff et al. 2008" ~ "Western Ontario", 
    Study == "Herbert-Read et al. 2019" ~ "Cambridge",  # Lund too
    Study == "Stephenson et al. 2015" ~ "Cardiff",
    Study == "Edenbrow et al. 2017" ~ "Exeter",
    Study == "Evans et al. 2003" ~ "Padova",
    Study == "Piyapong et al. 2011" ~ "Mahasarakham",
    Study == "Magurran and Seghers 1994" ~ "Oxford",
    Study == "Martin and Johnsen 2017" ~ "Duke",
    Study == "Gotanda et al. 2013" ~ "McGill",  ## 
    Study == "Hain et al. 2016" ~ "Western Ontario",
    Study == "de Lira et al. inprep" ~ "McGill",
    Study == "Eastya et al. 2011" ~ "McGill",
    Study == "Gordon et al. 2017b" ~ "Jyvaskyla",
    Study == "Fraser et al. 2015" ~ "Max Plank",
    Study == "Devigili et al. 2019" ~ "Stockholm",
    Study == "Perez-Jvostov et al. 2012" ~ "McGill",
    Study == "Dial et al. 2016" ~ "Brown",
    Study == "Egset et al. 2011" ~ "NTNU",
    Study == "El-Sabaawi et al. 2012" ~ "Cornell",
    Study == "Huizinga et al. 2009" ~ "Colorado State",
    Study == "Reznick et al. 2004" ~ "UC Riverside",
    Study == "Herbert-Read et al. 2017" ~ "Stockholm",
    Study == "Schwartz and Hendry 2007" ~ "McGill",
    Study == "Zandona et al. 2017" ~ "Drexel",
    Study == "Millar and Hendry 2011" ~ "McGill",
    Study == "Ioannou et al. 2017" ~ "Bristol",
    Study == "Reznick and Endler 1982" ~ "Pennsylvania",
    Study == "Reddon et al. 2018" ~ "McGill",
    Study == "Fischer et al. 2013" ~ "Colorado State",
    Study == "Elgee et al. 2010" ~ "Windsor",
    Study == "Gordon et al. 2012" ~ "UC Riverside",
    Study == "Evans and Magurran 1999" ~ "St Andrews",
    Study == "Dial et al. 2017" ~ "Brown",
    Study == "Bohr Brask et al. 2019" ~ "Exeter",
    Study == "Fischer et al. 2016" ~ "Colorado State",
    Study == "Valvo et al. 2019" ~ "Florida State",
    Study == "Reznick and Bryant 2007" ~ "UC Riverside",
    Study == "Magurran and Seghers" ~ "Oxford",
    Study == "Auer et al. 2018" ~ "Glasgow",
    Study == "Croft et al. 2009" ~ "Exeter",  ## also Bangor
    Study == "Zandona et al. 2015" ~ "Drexel"))  ## diff from current

spreadsheet.data <- 
  spreadsheet.data %>% 
  mutate(institution_NS = case_when(
    #Study == "O'Steen et al. 2002" ~ "UC Irvine",  
    #Study == "Weese et al. 2010" ~ "Maine",  
    Study == "Burns et al. 2009" ~ "Toronto",  
    Study == "Reznick et al. 2005" ~ "UC Riverside",
    #Study == "Evans et al. 2011" ~ "Western Australia",
    #Study == "Kam et al. 2014" ~ "Simon Fraser",  # This is a typo in spreadsheet, actually Sandkam
    Study == "Burns and Rodd 2008" ~ "Toronto",  
    #Study == "Harris et al. 2010" ~ "Lund",
    Study == "Neff et al. 2008" ~ "Western Ontario", 
    #Study == "Herbert-Read et al. 2019" ~ "Cambridge",  # Lund too
    #Study == "Stephenson et al. 2015" ~ "Cardiff",
    Study == "Edenbrow et al. 2017" ~ "Exeter",
    #Study == "Evans et al. 2003" ~ "Padova",
    #Study == "Piyapong et al. 2011" ~ "Mahasarakham",
    Study == "Magurran and Seghers 1994" ~ "Oxford",
    #Study == "Martin and Johnsen 2017" ~ "Duke",
    Study == "Gotanda et al. 2013" ~ "McGill",  ## 
    Study == "Hain et al. 2016" ~ "Western Ontario",
    Study == "de Lira et al. inprep" ~ "McGill",
    Study == "Eastya et al. 2011" ~ "McGill",
    #Study == "Gordon et al. 2017b" ~ "Jyvaskyla",
    #Study == "Fraser et al. 2015" ~ "Max Plank",
    Study == "Devigili et al. 2019" ~ "Stockholm",
    Study == "Perez-Jvostov et al. 2012" ~ "McGill",
    Study == "Dial et al. 2016" ~ "Brown",
    #Study == "Egset et al. 2011" ~ "NTNU",
    #Study == "El-Sabaawi et al. 2012" ~ "Cornell",
    Study == "Huizinga et al. 2009" ~ "Colorado State",
    Study == "Reznick et al. 2004" ~ "UC Riverside",
    Study == "Herbert-Read et al. 2017" ~ "Stockholm",
    Study == "Schwartz and Hendry 2007" ~ "McGill",
    Study == "Zandona et al. 2017" ~ "Drexel",
    Study == "Millar and Hendry 2011" ~ "McGill",
    #Study == "Ioannou et al. 2017" ~ "Bristol",
    #Study == "Reznick and Endler 1982" ~ "Pennsylvania",
    Study == "Reddon et al. 2018" ~ "McGill",
    Study == "Fischer et al. 2013" ~ "Colorado State",
    #Study == "Elgee et al. 2010" ~ "Windsor",
    Study == "Gordon et al. 2012" ~ "UC Riverside",
    #Study == "Evans and Magurran 1999" ~ "St Andrews",
    Study == "Dial et al. 2017" ~ "Brown",
    Study == "Bohr Brask et al. 2019" ~ "Exeter",
    Study == "Fischer et al. 2016" ~ "Colorado State",
    #Study == "Valvo et al. 2019" ~ "Florida State",
    Study == "Reznick and Bryant 2007" ~ "UC Riverside",
    Study == "Magurran and Seghers" ~ "Oxford",
    #Study == "Auer et al. 2018" ~ "Glasgow",
    Study == "Croft et al. 2009" ~ "Exeter",  ## also Bangor
    Study == "Zandona et al. 2015" ~ "Drexel"))  ## diff from current

## QUESTION ECOLOGY (NORTH VS SOUTH)

## combine R2 and spreadsheet data
#prep for when we bind them together, this variable will go in the model to indicate the type of R2
R2.data.among$method <- "all"
R2.data.south$method <- "south"
R2.data.intro$method <- "intro"
R2.data.caroni$method <- "caroni"
R2.data.among.drainage$method <- "both.drainages"

#get relevant data with one entry for each Trait
##I was inclusive with columns, many probably aren't needed so you can cut them out
data.all <- inner_join(spreadsheet.data, R2.data.among, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:54)%>% #I selected only the columns that have information that applies at the trait level
  distinct(TraitID, .keep_all = TRUE) #removes replicated TraitIDs and retains the columns

#repeat for south only R2 and others...
data.south <- inner_join(spreadsheet.data, R2.data.south, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:54)%>% 
  distinct(TraitID, .keep_all = TRUE) 

data.intro<- inner_join(spreadsheet.data, R2.data.intro, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:54)%>% 
  distinct(TraitID, .keep_all = TRUE) 

data.caroni<- inner_join(spreadsheet.data, R2.data.caroni, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:54)%>% 
  distinct(TraitID, .keep_all = TRUE) 

data.among.drainage<- inner_join(spreadsheet.data, R2.data.among.drainage, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:54)%>% 
  distinct(TraitID, .keep_all = TRUE)

#then we can bind them together as we wish! First ecology...
data.for.ecology.models<-rbind(data.all,data.south) %>% 
  arrange(TraitID) %>% #puts them in a nice order
  group_by(TraitID) %>% #groups them for the count
  filter(n() > 1) #filters only trait IDs that have more than 1 entry within the group (n = 204 traits)

ecology.data.males <- data.for.ecology.models %>% 
  filter(Sex == "M" & Kingsolver_traits !="Other")
ecology.data.females <- data.for.ecology.models %>% 
  filter(Sex == "F" & Kingsolver_traits !="Other")

#evolutionary history question
data.for.evolhist.models<-rbind(data.caroni,data.among.drainage) %>% 
  arrange(TraitID) %>% #puts them in a nice order
  group_by(TraitID) %>% #groups them for the count
  filter(n() > 1) #filters only trait IDs that have more than 1 entry within the group (n = 204 traits)

evolhist.data.males <- data.for.evolhist.models %>% 
  filter(Sex == "M" & Kingsolver_traits !="Other")
evolhist.data.females <- data.for.evolhist.models %>% 
  filter(Sex == "F" & Kingsolver_traits !="Other")

# Intro, I'm pretty sure these are the right csvs for this question
data.for.intro.models<-rbind(data.all,data.intro) %>% 
  arrange(TraitID) %>% #puts them in a nice order
  group_by(TraitID) %>% #groups them for the count
  filter(n() > 1) #filters only trait IDs that have more than 1 entry within the group (n = 204 traits)

time.data.males <- data.for.intro.models %>% 
  filter(Sex == "M" & Kingsolver_traits !="Other")
time.data.females <- data.for.intro.models %>% 
  filter(Sex == "F" & Kingsolver_traits !="Other")

# as a note and fail safe I would run the grouping with TraitID on the sex specific dfs to make sure theres two traitID entries for each

## Overall models

## remove 'both'
### (Because when calculated by hand, I would do M/F/Both with same data)
data.all <- data.all %>% filter(Sex %in% c("M", "F"))  

## remove other (because overall trait type model will not converge with other)
data.all <- data.all %>% filter(!Kingsolver_traits == "Other")

## creating a dataset w no colour (to compare w those w colour)
data.all.no.colour <- data.all %>% filter(!Kingsolver_traits == "Colour")

## basic histograms
data.all %>% 
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

## overall traits model

## this WILL NOT RUN with other included
(all.model.traits <- glmer(R.2 ~ Kingsolver_traits +  (1|StudyID), data = data.all, family = binomial)) %>% summary()

## overall sex model
(all.model.sex <- glmer(R.2 ~ Sex + (1|StudyID), data = data.all, family = binomial)) %>% summary()

### Removed colour 
(all.model.sex <- glmer(R.2 ~ Sex + (1|StudyID), data = data.all.no.colour, family = binomial)) %>% summary()

## trait type plots
## here is a general plot with colour
data.all %>% 
  ggplot(aes(x = Kingsolver_traits, y = R.2, fill = Kingsolver_traits)) + 
  geom_jitter(aes(color = Kingsolver_traits), width = .1) +
  geom_boxplot(size = 1, aes(x = Kingsolver_traits, y = R.2) ,
               alpha = 0.3, width = 0.3) +
  scale_colour_manual(values =  c("#FFB93C", "#457111","#15899A", "#BD8DC3", "#D07D7D", "dark gray")) +
  scale_fill_manual(values =  c("#FFB93C", "#457111","#15899A", "#BD8DC3", "#D07D7D", "dark gray")) +
  theme_classic() 

data.all.no.colour %>% 
  ggplot(aes(x = Kingsolver_traits, y = R.2, fill = Kingsolver_traits)) + 
  geom_jitter(aes(color = Kingsolver_traits), width = .1) +
  geom_boxplot(size = 1, aes(x = Kingsolver_traits, y = R.2) ,
               alpha = 0.3, width = 0.3) +
  scale_colour_manual(values =  c("#FFB93C", "#457111","#15899A", "#BD8DC3", "#D07D7D", "dark gray")) +
  scale_fill_manual(values =  c("#FFB93C", "#457111","#15899A", "#BD8DC3", "#D07D7D", "dark gray")) +
  theme_classic() 

## sex plots
data.all %>% 
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

# Models ----

## 1. Ecology models ----
#### Is there a difference between the slopes? 

### Remove 'Both' sex category because duplicates
data.for.ecology.models <- data.for.ecology.models %>% filter(Sex %in% c("M", "F"))  

### Remove other (to be consistent with above, but idk)
data.for.ecology.models %>% filter(!Kingsolver_traits == "Other")

### Make dataframes to remove colour
ecology.data.no.colour <- data.for.ecology.models %>% filter(!Kingsolver_traits == "Colour")

### Model with everything
(ecology.full <- glmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.ecology.models, family = binomial)) %>% summary()

### Model without colour
(ecology.no.colour <- glmer(R.2 ~ method*Sex + (1|StudyID), data = ecology.data.no.colour, family = binomial)) %>% summary()

### effect size plot with all models
plot_models(ecology.full, ecology.no.colour, vline.color = "grey")

### this is for the means that I report in the text
ecology.data.no.colour %>% filter(method == "all") %>% summary()
ecology.data.no.colour %>% filter(method == "south") %>% summary()

### Plots

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


## 2. Intro models ----
#### Is there a difference between studies with only natural vs w intro

### Remove 'Both' sex category because duplicates
data.for.intro.models <- data.for.intro.models %>% filter(Sex %in% c("M", "F"))  

### Remove 'other' to be consistent
data.for.intro.models <- data.for.intro.models %>% filter(!Kingsolver_traits == "Other")

### Make dataframes to remove colour
intro.data.no.colour <- data.for.intro.models %>% filter(!Kingsolver_traits == "Colour")

### Model with everything 
(intro.full <- glmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.intro.models, family = binomial)) %>% summary()

### Model without colour
(intro.no.colour <- glmer(R.2 ~ method*Sex + (1|StudyID), data = intro.data.no.colour, family = binomial)) %>% summary()

### effect size plot with all models
plot_models(intro.full, intro.no.colour, vline.color = "grey")

intro.data.no.colour %>% 
  ggplot(aes(x = method, y = R.2, fill = Sex)) +
  geom_boxplot(size = 1.5) +
  geom_jitter(size = 3, alpha = 0.1, width = 0.2) +
  scale_fill_manual(values = c("gray48", "gray90")) +
  facet_wrap(~Sex) +
  theme_bw()

### this is for the means that I report in the text
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

## 3. Evolutionary history models ----
#### Is there a difference when only w pops in the caroni vs also in the Oropuche? 

### These are all singular fits

### Remove 'Both' sex category because duplicates
data.for.evolhist.models <- data.for.evolhist.models %>% filter(Sex %in% c("M", "F"))  

### Remove other, being consistent
data.for.evolhist.models <- data.for.evolhist.models %>% filter(!Kingsolver_traits == "Other")

### Make dataframes to remove colour
evolhist.data.no.colour <- data.for.evolhist.models %>% filter(!Kingsolver_traits == "Colour")

### Model with everything 
(evolhist.full <- glmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.evolhist.models, family = binomial)) %>% summary()

### Model without colour
(evolhist.no.colour <- glmer(R.2 ~ method*Sex + (1|StudyID), data = evolhist.data.no.colour, family = binomial)) %>% summary()

### effect size plot with all models
plot_models(evolhist.full, evolhist.no.colour, vline.color = "grey")

### These are means I report in text
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


# Troubleshooting ----
#### Evolhist are singular, so here we are comparing our models that are not singular to glms/lmer, to see if it changes anything.
#### Note that I can't plot lmer w an sjplot (cause lmer is plotted as an Estimate, whereas the glm/glmer are odds ratios). I think that the glm/glmer are more consistent anyway (?) but I am only basing this on agreement re significance.

## Compare intro
### This one is not singular (so is an example)
summary(intro.no.colour)
(intro.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = intro.data.no.colour)) %>% summary()
Anova(intro.full.lmer, type = 2)
(intro.full.glm <- glm(R.2 ~ method*Sex, data = intro.data.no.colour, family = binomial)) %>% summary()

### Plot models
plot_models(intro.no.colour, intro.full.glm)  # I can't tell

## Compare ecology
summary(ecology.no.colour)
(ecology.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = ecology.data.no.colour)) %>% summary()
Anova(ecology.full.lmer, type = 3)
(ecology.full.glm <- glm(R.2 ~ method*Sex, data = ecology.data.no.colour, family = binomial)) %>% summary()

### Plot model to compare effect sizes
plot_models(ecology.no.colour, ecology.full.glm)

### This is singular (need to decide with this one)
summary(evolhist.no.colour)
(evolhist.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = evolhist.data.no.colour)) %>% summary()
Anova(evolhist.full.lmer, type = 3)
(evolhist.full.glm <- glm(R.2 ~ method*Sex, data = evolhist.data.no.colour, family = binomial)) %>% summary()

plot_models(evolhist.no.colour, evolhist.full.glm)

### Plot models (is there a better way to compare effect sizes?)




# !!!
# !!!
# Everything below here is old!!!!
## models

## ECOLOGY QUESTION (ACROSS SLOPES VS WITHIN SOUTH SLOPE)

## glmer first
(ecology.model.males.glmer <- glmer(R.2 ~ method + (1|institution), family = binomial,
                                    data = ecology.data.males)) %>% summary()

(ecology.model.females.glmer <- glmer(R.2 ~ method + (1|institution), family = binomial,
                                      data = ecology.data.females)) %>% summary()

(ecology.model.males.glmer <- glmer(R.2 ~ method + (1|institution_NS), family = binomial,
                                    data = ecology.data.males)) %>% summary()

### this is the only one without the warning
(ecology.model.females.glmer <- glmer(R.2 ~ method + (1|institution_NS), family = binomial,
                                      data = ecology.data.females)) %>% summary()

(ecology.model.males.glmer <- glmer(R.2 ~ method + (1|first_author), family = binomial,
                                    data = ecology.data.males)) %>% summary()

(ecology.model.females.glmer <- glmer(R.2 ~ method + (1|first_author), family = binomial,
                                      data = ecology.data.females)) %>% summary()

(ecology.model.males.glmer <- glmer(R.2 ~ method + (1|first_author_NS), family = binomial,
                                    data = ecology.data.males)) %>% summary()

(ecology.model.females.glmer <- glmer(R.2 ~ method + (1|first_author_NS), family = binomial,
                                      data = ecology.data.females)) %>% summary()

## glm second

(ecology.model.males.glm <- glm(R.2 ~ method + institution, family = binomial,
                                data = ecology.data.males)) %>% summary()

(ecology.model.females.glm <- glm(R.2 ~ method + institution, family = binomial,
                                  data = ecology.data.females)) %>% summary()

(ecology.model.males.glm <- glm(R.2 ~ method + first_author, family = binomial,
                                data = ecology.data.males)) %>% summary()

(ecology.model.females.glm <- glm(R.2 ~ method + first_author, family = binomial,
                                  data = ecology.data.females)) %>% summary()

(ecology.model.males.glm <- glm(R.2 ~ method + StudyID, family = binomial,
                                data = ecology.data.males)) %>% summary()

(ecology.model.females.glm <- glm(R.2 ~ method + StudyID, family = binomial,
                                  data = ecology.data.females)) %>% summary()


## TIME FRAME QUESTION (WITH INTORS VS ONLY NATURAL)
## glmer first
(time.model.males.glmer <- glmer(R.2 ~ method + (1|institution), family = binomial,
                                 data = time.data.males)) %>% summary()

### no warning
(time.model.females.glmer <- glmer(R.2 ~ method + (1|institution), family = binomial,
                                   data = time.data.females)) %>% summary()

### no warning
(time.model.males.glmer <- glmer(R.2 ~ method + (1|institution_NS), family = binomial,
                                 data = time.data.males)) %>% summary()

### no warning
(time.model.females.glmer <- glmer(R.2 ~ method + (1|institution_NS), family = binomial,
                                   data = time.data.females)) %>% summary()

(time.model.males.glmer <- glmer(R.2 ~ method + (1|first_author), family = binomial,
                                 data = time.data.males)) %>% summary()

### no warning
(time.model.females.glmer <- glmer(R.2 ~ method + (1|first_author), family = binomial,
                                   data = time.data.females)) %>% summary()

### no warning
(time.model.males.glmer <- glmer(R.2 ~ method + (1|first_author_NS), family = binomial,
                                 data = time.data.males)) %>% summary()

### no warning
(time.model.females.glmer <- glmer(R.2 ~ method + (1|first_author_NS), family = binomial,
                                   data = time.data.females)) %>% summary()

## glm second


(time.model.males.glm <- glm(R.2 ~ method + institution, family = binomial,
                             data = time.data.males)) %>% summary()

(time.model.females.glm <- glm(R.2 ~ method + institution, family = binomial,
                               data = time.data.females)) %>% summary()

(time.model.males.glm <- glm(R.2 ~ method + first_author, family = binomial,
                             data = time.data.males)) %>% summary()

(time.model.females.glm <- glm(R.2 ~ method + first_author, family = binomial,
                               data = time.data.females)) %>% summary()

# QUESTION EVOLUTIONARY HISTORY (CARONI VS OROPUCHE)

## glmer first
(evolhist.model.males.glmer <- glmer(R.2 ~ method + (1|institution), family = binomial,
                                     data = evolhist.data.males)) %>% summary()

(ecology.model.females.glmer <- glmer(R.2 ~ method + (1|institution), family = binomial,
                                      data = evolhist.data.females)) %>% summary()

(evolhist.model.males.glmer <- glmer(R.2 ~ method + (1|institution_NS), family = binomial,
                                     data = evolhist.data.males)) %>% summary()

(ecology.model.females.glmer <- glmer(R.2 ~ method + (1|institution_NS), family = binomial,
                                      data = evolhist.data.females)) %>% summary()

(evolhist.model.males.glmer <- glmer(R.2 ~ method + (1|first_author), family = binomial,
                                     data = evolhist.data.males)) %>% summary()

(ecology.model.females.glmer <- glmer(R.2 ~ method + (1|first_author), family = binomial,
                                      data = evolhist.data.females)) %>% summary()

(evolhist.model.males.glmer <- glmer(R.2 ~ method + (1|first_author_NS), family = binomial,
                                     data = evolhist.data.males)) %>% summary()

(ecology.model.females.glmer <- glmer(R.2 ~ method + (1|first_author_NS), family = binomial,
                                      data = evolhist.data.females)) %>% summary()

## glm second


(evolhist.model.males.glm <- glm(R.2 ~ method + institution, family = binomial,
                                 data = evolhist.data.males)) %>% summary()

(evolhist.model.females.glm <- glm(R.2 ~ method + institution, family = binomial,
                                   data = evolhist.data.females)) %>% summary()

(evolhist.model.males.glm <- glm(R.2 ~ method + first_author, family = binomial,
                                 data = evolhist.data.males)) %>% summary()

(evolhist.model.females.glm <- glm(R.2 ~ method + first_author, family = binomial,
                                   data = evolhist.data.females)) %>% summary()


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
