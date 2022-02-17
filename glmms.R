# GLMMs for parallel_guppies! 
# 2022-02-09 cleaned up

# Libraries ---- 
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
library(nlme)


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

## remove 'both'
(data.all <- data.all %>% filter(Sex %in% c("M", "F"))  ) %>% summary()

# Overall models (traits, sex, rearing) ----

## remove other (because model below will not converge with other)

## Trait type model (in paper) ----
data.all.traits <- data.all %>% filter(!Kingsolver_traits == "Other")
(all.model.traits <- glmer(R.2 ~ Kingsolver_traits +  (1|StudyID), 
                           data = data.all.traits, family = binomial)) %>% summary()
car::Anova(all.model.traits, type = "II")


## sex with colour (in paper) ----
## with colour
(all.model.sex <- glmer(R.2 ~ Sex + (1|StudyID), data = data.all, family = binomial)) %>% summary()
car::Anova(all.model.sex, type = "II")

## sex without colour  (in paper) ----
(all.model.sex.no.colour <- glmer(R.2 ~ Sex + (1|StudyID), data = data.all.no.colour, family = binomial)) %>% summary()
car::Anova(all.model.sex, type = "II")

## Rearing enviro model (in paper) ----
data.all.rear <- data.all %>% filter(StudyType %in% c("Common Garden (F2)", "Wildcaught"))  # won't run w CG F1 (not a lot anyway)
(all.model.rearing <- glmer(R.2 ~ StudyType +  (1|StudyID), data = data.all.rear, family = binomial)) %>% summary()
car::Anova(all.model.rearing, type = "II")


# multivariate models (traits, sex, rearing) ----

## sex and traits (in paper) ----
(sex.and.triats <- glmer(R.2 ~ Kingsolver_traits + Sex + (1|StudyID), data = data.all, family = binomial)) %>% summary()
Anova(sex.and.triats, type = "II")

## sex and rear (in paper) ----
(sex.and.rear <- glmer(R.2 ~ StudyType + Sex + (1|StudyID), data = data.all.rear, family = binomial)) %>% summary()
Anova(sex.and.rear, type = 2)

# Determinants models ----
## Ecology models ----

### fix structure
data.for.ecology.models$method <- as.factor(data.for.ecology.models$method)

### Remove 'Both' sex category because duplicates
data.for.ecology.models <- data.for.ecology.models %>% filter(Sex %in% c("M", "F"))  

## Ecology model (in paper)
(ecology.full <- glmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.ecology.models, family = binomial)) %>% summary()
car::Anova(ecology.full, type = "II")

## remove the interaction (in paper)
(ecology.full <- glmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.ecology.models, family = binomial)) %>% summary()
car::Anova(ecology.full, type = "II")

## Intro models ----
# (Note that I deleted the old intro file - 2022-02-09 - this is only intro broad)
## Fix structure
data.for.intro.models.broad$method <- as.factor(data.for.intro.models.broad$method)

### Remove 'Both' sex category because duplicates
data.for.intro.models.broad <- data.for.intro.models.broad %>% filter(Sex %in% c("M", "F"))  

## Intro model (in paper)
(intro.full.broad <- glmer(R.2 ~ method + (1|StudyID), data = data.for.intro.models.broad, family = binomial)) %>% summary()
car::Anova(intro.full.broad, type = "II")


## Evolutionary history models ----
# (These are all singular fits)

## fix structure
data.for.evolhist.models$method <- as.factor(data.for.evolhist.models$method)

## Remove 'Both' sex category because duplicates
data.for.evolhist.models <- data.for.evolhist.models %>% filter(Sex %in% c("M", "F"))  

## Evolutionary history model (in paper)
(evolhist.full <- glmer(R.2 ~ method + Sex + (1|StudyID), data = data.for.evolhist.models, family = binomial)) %>% summary()
car::Anova(evolhist.full, type = "II")

(evolhist.glm <- glm(R.2 ~ method + StudyID, data = data.for.evolhist.models, family = binomial)) %>% summary()
car::Anova(evolhist.glm, type = "II")


### Troubleshooting Evolutionary history----
### Evolhist are all singular, so here we are comparing our models that are not singular to glms/lmer, to see if it changes anything.

### This is singular (need to decide with this one)

summary(evolhist.full) # GLMM above

# try w lmer
(evolhist.full.lmer <- lmer(R.2 ~ method*Sex + (1|StudyID), data = data.for.evolhist.models)) %>% summary()
Anova(evolhist.full.lmer, type = 3)

# tri w glm
(evolhist.full.glm <- glm(R.2 ~ method*Sex, data = data.for.evolhist.models, family = binomial)) %>% summary()

car::Anova(evolhist.full.lmer, type = "II")
car::Anova(evolhist.full.glm, type = "II")

# Other means etc reported in text ----

# Traits
(Colourtrait <- data.all %>% filter(Kingsolver_traits == "Colour")) %>% summary()
(LHtrait <- data.all %>% filter(Kingsolver_traits == "Other_life_history")) %>% summary()
(Sizetrait <- data.all %>% filter(Kingsolver_traits == "Size")) %>% summary()
(Morphtrait <- data.all %>% filter(Kingsolver_traits == "Other_morphology")) %>% summary()
(Othertrait <- data.all %>% filter(Kingsolver_traits == "Other")) %>% summary()
(Phystrait <- data.all %>% filter(Kingsolver_traits == "Physiology")) %>% summary()
(Behavtrait <- data.all %>% filter(Kingsolver_traits == "Behaviour")) %>% summary()

# Sex
(Malewithcolour <- data.all %>% filter(Sex == "M")) %>% summary()
(Malenocolour <- data.all.no.colour %>% filter(Sex == "M")) %>% summary()
(Femalesex <- data.all %>% filter(Sex == "F")) %>% summary()

# Rearing enviro
(cg <- data.all %>% filter(StudyType == "Common Garden (F2)")) %>% summary()
(wc <- data.all %>% filter(StudyType == "Wildcaught")) %>% summary()

## Ecology
(one.slope <- data.for.ecology.models %>% filter(method == "south")) %>% summary()
(both.slopes <- data.for.ecology.models %>% filter(method == "all")) %>% summary()


## Intro
(only_natural_broad <- data.for.intro.models.broad %>% filter(method == "only_natural_broad")) %>% summary()
(natural_and_intro_broad <- data.for.intro.models.broad %>% filter(method == "all")) %>% summary()

## Evolhist
(one_drainage <- data.for.evolhist.models %>% filter(method == "caroni")) %>% summary()
(both_drainages <- data.for.evolhist.models %>% filter(method == "both.drainages")) %>% summary()


# figures in manuscript ----

## figure 1 is map 

## figure 2 ----
fig2data <- read.csv("testPG_figure.csv", fileEncoding="UTF-8-BOM")
fig2data$population <- as.factor(fig2data$population)

fig2data$population <- factor(fig2data$population, levels = (c("HP1", "LP1", "HP2", "LP2")))

levels(fig2data$population)

highparal <- fig2data %>% filter(facet == 'a') %>% 
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
        plot.title = element_text(size = 16))


lowparal <- fig2data %>% filter(facet == 'b') %>% 
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
        plot.title = element_text(size = 16))

ggarrange(highparal, lowparal, common.legend = TRUE, legend = "bottom")

(highparal <- highparal + annotate("segment", x = 1, xend = 2, y = 12.01, yend = 13.5))
(highparal <- highparal + annotate("segment", x = 3, xend = 4, y = 12.1, yend = 13.76))
highparal

(lowparal <- lowparal + annotate("segment", x = 1, xend = 2, y = 8.74, yend = 9.85))
(lowparal <- lowparal + annotate("segment", x = 3, xend = 4, y = 8.2, yend = 7.0))
lowparal

fig2 <- ggarrange(highparal, lowparal, common.legend = TRUE, legend = "bottom")

#tiff("pg.fig2.tiff", res = 600, units = "in", height = 5, width = 7)
fig2
#dev.off()

## figure 3 ----

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

m_hist_w_colour <- data.all %>% filter(Sex %in% c("M", "F"))
m_hist_w_colour <- 
  m_hist_w_colour %>% mutate(renamedSex =
                           case_when(Sex == "M" ~ "Male (with colour)",
                              Sex == "F" ~ "Female"))

m_hist <- 
  m_hist_w_colour %>% 
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
  facet_wrap(. ~ renamedSex) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

m_hist_no_colour <- data.all %>% filter(Sex %in% c("M", "F") &
                                          !Kingsolver_traits == "Colour")
  
m_hist_no_colour <- 
  m_hist_no_colour %>% mutate(renamedSex =
                               case_when(Sex == "M" ~ "Male (without colour)",
                                         Sex == "F" ~ "Female"))


m_hist_no_colour <- 
  m_hist_no_colour %>% 
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
  facet_wrap(. ~ renamedSex) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

f_hist_w_other <-
  m_hist_w_colour %>% 
  filter(Sex == "F") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 #fill="#E6E6E6", 
                 colour = "black", fill = "#E6E6E6", size = 1) +
  #xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  # ggtitle("Females (n = 172)") +
  facet_wrap(. ~ renamedSex) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

top.fig.3 <- cowplot::plot_grid("", overall_hist + theme(axis.title = element_blank()),
                                wildcaught_hist + theme(axis.title = element_blank()),
                                cg_hist + theme(axis.title = element_blank()), 
                                nrow = 1, labels = c("", "A", "B", "C"), 
                                rel_widths = c(0.15, 1, 1, 1))

bottom.fig.3 <- cowplot::plot_grid("", m_hist + theme(axis.title = element_blank()),
                                   m_hist_no_colour + theme(axis.title.y = element_blank()),
                                   f_hist_w_other + theme(axis.title = element_blank()),
                                   
                                   #f_hist_no_other + theme(axis.title.y = element_blank()), 
                                   nrow = 1,
                                   labels = c("", "D", "E", "F", "G"),
                                   align = "h", axis = "bt",
                                   rel_widths = c(0.15, 1, 1,1))

fig3 <- cowplot::plot_grid(top.fig.3, bottom.fig.3, nrow = 2, rel_heights = c(1,1))
fig3 <- figure3 + cowplot::draw_label("Frequency", x=  0, y=0.5, vjust= 1, angle=90, size = 24)

#tiff("pg.figure3.tiff", res= 600, units = "in", height = 6, width = 8)
fig3
#dev.off()

## figure 4 ----
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
fig4 <- cowplot::plot_grid(toptraits,bottomtraits, nrow = 2, align = "bt", rel_heights = c(0.9, 1))
fig4 <- traitsplot + cowplot::draw_label("Frequency", x=  0, y=0.5, vjust= 1, angle=90, size = 24)

#tiff("fig4.pg.tiff", res = 600, width = 8, height = 5, units = "in")
fig4
#dev.off()

## figure 5 ----

### intro
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

### ecology

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

### evolutionary history

data.evo <- data.for.evolhist.models
levels(data.evo$method) <- c("Between both drainages",
                             "Only the Caroni drainage")
(evolhist_hist_a <-
    data.evo %>% 
    filter(method == "Between both drainages") %>% 
    ggplot(aes(x = R.2)) +
    geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                   fill="#E6E6E6", colour = "black", size = 1) +
        facet_wrap(~method) +
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
#tiff("pg.fig5.tiff", res= 600, units = "in", height = 6, width =8)
fig5
dev.off()

## figure 6 ----

# (These dataframes are all made above)

listforfig6 <- list("Male with colour (274)" = Malewithcolour$R.2, "Male without colour (165)" = Malenocolour$R.2, 
           "Female (172)" = Femalesex$R.2, 
           "Colour (109)" = Colourtrait$R.2, "Life history (48)" = LHtrait$R.2, "Size (47)" = Sizetrait$R.2, 
           "Morphology (42)" = Morphtrait$R.2, "Other (3)" = Othertrait$R.2, "Physiology (64)" = Phystrait$R.2,
           "Common Garden (F2) (70)" = cg$R.2, "Wild caught (373)" = wc$R.2, "Within one slope (176)" = one.slope$R.2, "Between both slopes (176)" = both.slopes$R.2, 
           "Only natural (183)" = only_natural_broad$R.2, "Natural and introduced (183)" = natural_and_intro_broad$R.2, 
           "Within one drainage (258)" = one_drainage$R.2, "Between both drainages (258)" = both_drainages$R.2)

mean <- (as.data.frame(sapply(listforfig6, mean, na.rm = TRUE)))
mean <- cbind(Factor = rownames(mean), mean)
rownames(mean) <- NULL
colnames(mean)[2] <- "mean"

standev <- (as.data.frame(sapply(listforfig6, sd, na.rm = TRUE)))
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

(figure6<-
    larger_table %>% ggplot(x = Factor, aes(reorder(Factor, order), y = mean)) +
    geom_point(size = 2) +
    geom_linerange(aes(x = Factor, ymin = minsd, ymax = maxsd), size = 0.5) +
    
    annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf, fill = "yellow", alpha = 0.2) +
    annotate("rect", xmin = 3.5, xmax = 9.5, ymin = -Inf, ymax = Inf, fill = "orange", alpha = .4) +
    annotate("rect", xmin = 9.5, xmax = 11.5, ymin = -Inf, ymax = Inf, fill = "palegreen", alpha = .4) +
    annotate("rect", xmin = 11.5, xmax = 13.5, ymin = -Inf, ymax = Inf, fill = "#0404B7", alpha = .3) +
    annotate("rect", xmin = 13.5, xmax = 15.5, ymin = -Inf, ymax = Inf, fill = "purple", alpha = .4) +
    annotate("rect", xmin = 15.5, xmax = 17.5, ymin = -Inf, ymax = Inf, fill = "hot pink", alpha = .3) +
    
    
    labs(x = "", 
         y = expression(paste(R^2))) +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 24)))

tiff("pg.figure6.tiff", res = 600, units = "in", height = 6, width = 8)
figure6
dev.off()


# Supplementary material ----
## Supplemental table 1 ----

suptable.R2 <- data.all %>% group_by(Study)
suptable.R2 <- mutate(suptable.R2, min(R.2))
suptable.R2 <- mutate(suptable.R2, max(R.2))

suptable.people <- data.frame(data.all %>% group_by(Study, Poptype_broad) %>% tally())

binded.suptable <- rbind(suptable.R2, suptable.people)

No.traits <- data.all %>% group_by(Study, Trait) %>% tally() %>% 
  transmute(No.traits = sum(n)) %>% 
  distinct(No.traits, .keep_all = FALSE) 

Slope.names <- spreadsheet.data %>% group_by(Study, Slope) %>% count()
Slope.names <- data.frame(aggregate(Slope ~ Study, data = Slope.names, paste, collapse = ","))

traits.slopes <- left_join(No.traits, Slope.names)

Drainage.names <- spreadsheet.data %>% group_by(Study, Drainage) %>% count()
Drainage.names <- data.frame(aggregate(Drainage ~ Study, data = Drainage.names, paste, collapse = ","))

traits.slopes.drainages <- left_join(traits.slopes, Drainage.names)

Pop.type <- spreadsheet.data %>% group_by(Study, Poptype_broad) %>% count()
Pop.type <- data.frame(aggregate(Poptype_broad ~ Study, data = Pop.names, paste, collapse = ","))

traits.slopes.drinages.pops <- left_join(traits.slopes.drainages, Pop.type)


Populations.names <- spreadsheet.data %>% group_by(Study, Poptype_broad, Population) %>% 
  filter(Poptype_broad == "Introduction") %>% count()
Populations.names <- data.frame(aggregate(Population ~ Study, data = Populations.names, paste, collapse = ","))

final.sup.table <- left_join(traits.slopes.drinages.pops, Populations.names)

View(final.sup.table)

## Supplemental figure 1 ----

(sup1<-
   data.all %>% 
   filter(!Kingsolver_traits %in% c("Colour", "Other")) %>% 
   ggplot(aes(y = R.2)) +
   geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                  fill="#E6E6E6", colour = "black", size = 1) +
   xlab(expression(paste(R^2))) +
   ylab("Frequency (%)") +
   theme_bw() +
   facet_grid(Sex ~ Kingsolver_traits, scales = "free") +
   theme(
     strip.text = element_text(size = 13),
     axis.title = element_text(size = 24)))

behav_hist_f <-
  data.all %>% 
  filter(Kingsolver_traits == "Behaviour" &
           Sex == "F") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  xlim(0,1) +
  ylim(0,35)+
  #ggtitle("Wildcaught (n = 373)") +
  facet_wrap(. ~ Kingsolver_traits) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

behav_hist_m <-
  data.all %>% 
  filter(Kingsolver_traits == "Behaviour" &
           Sex == "M") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  xlim(0,1) +
  ylim(0,35)+
  #ggtitle("Wildcaught (n = 373)") +
  #facet_wrap(. ~ Kingsolver_traits) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

behav_sex_hist <- cowplot::plot_grid(behav_hist_f +theme(axis.title = element_blank(),
                                                         axis.text.x = element_blank()), 
                                     behav_hist_m + theme(axis.title = element_blank()), 
                                     ncol = 1)

behav_sex_hist

LH_hist_f <-
  data.all %>% 
  filter(Kingsolver_traits == "Other_life_history" &
           Sex == "F") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  xlim(0,1) +
  ylim(0,35)+
  #ggtitle("Wildcaught (n = 373)") +
  facet_wrap(. ~ Kingsolver_traits) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

LH_hist_m <-
  data.all %>% 
  filter(Kingsolver_traits == "Other_life_history" &
           Sex == "M") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  xlim(0,1) +
  ylim(0,35)+
  #ggtitle("Wildcaught (n = 373)") +
  #facet_wrap(. ~ Kingsolver_traits) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

LH_sex_hist <- cowplot::plot_grid(LH_hist_f +theme(axis.title = element_blank(),
                                                   axis.text.x = element_blank()), 
                                  LH_hist_m + theme(axis.title = element_blank()), 
                                  ncol = 1)

morph_hist_f <-
  data.all %>% 
  filter(Kingsolver_traits == "Other_morphology" &
           Sex == "F") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  xlim(0,1) +
  ylim(0,35)+
  #ggtitle("Wildcaught (n = 373)") +
  facet_wrap(. ~ Kingsolver_traits) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))


morph_hist_m <-
  data.all %>% 
  filter(Kingsolver_traits == "Other_morphology" &
           Sex == "M") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  xlim(0,1) +
  ylim(0,35)+
  #ggtitle("Wildcaught (n = 373)") +
  #facet_wrap(. ~ Kingsolver_traits) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

morph_sex_hist <- cowplot::plot_grid(morph_hist_f +theme(axis.title = element_blank(),
                                                         axis.text.x = element_blank()),
                                     morph_hist_m + theme(axis.title = element_blank()), 
                                     ncol = 1)

size_hist_f <-
  data.all %>% 
  filter(Kingsolver_traits == "Size" &
           Sex == "F") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  xlim(0,1) +
  ylim(0,35)+
  #ggtitle("Wildcaught (n = 373)") +
  facet_wrap(. ~ Kingsolver_traits) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

size_hist_m <-
  data.all %>% 
  filter(Kingsolver_traits == "Size" &
           Sex == "M") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  xlim(0,1) +
  ylim(0,35)+
  #ggtitle("Wildcaught (n = 373)") +
  #facet_wrap(. ~ Kingsolver_traits) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

size_sex_hist <- cowplot::plot_grid(size_hist_f +theme(axis.title = element_blank(),
                                                       axis.text.x = element_blank()),
                                    size_hist_m + theme(axis.title = element_blank()), 
                                    ncol = 1)

phys_hist_f <-
  data.all %>% 
  filter(Kingsolver_traits == "Physiology" &
           Sex == "F") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  xlim(0,1) +
  ylim(0,35)+
  #ggtitle("Wildcaught (n = 373)") +
  facet_wrap(. ~ Kingsolver_traits) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

phys_hist_m <-
  data.all %>% 
  filter(Kingsolver_traits == "Physiology" &
           Sex == "M") %>% 
  ggplot(aes(x = R.2)) +
  geom_histogram(mapping=aes(x=R.2, y=..count../sum(..count..)*100), bins=10, 
                 fill="#E6E6E6", colour = "black", size = 1) +
  xlab(expression(paste(R^2))) +
  ylab("Frequency") +
  xlim(0,1) +
  ylim(0,35)+
  #ggtitle("Wildcaught (n = 373)") +
  #facet_wrap(. ~ Kingsolver_traits) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 24))

phys_sex_hist <- cowplot::plot_grid(phys_hist_f +theme(axis.title = element_blank(),
                                                       axis.text.x = element_blank()),
                                    phys_hist_m + theme(axis.title = element_blank()), 
                                    ncol = 1)


sex_facet_hist <- cowplot::plot_grid("",behav_sex_hist, LH_sex_hist, morph_sex_hist, 
                                     phys_sex_hist, size_sex_hist, nrow = 1,
                                     rel_widths = c(0.20,1,1,1,1,1))

sex_facet_hist <- sex_facet_hist + cowplot::draw_label("Frequency (%)", x=  0, y=0.5, vjust= 1, angle=90, size = 24)

#tiff("figs1.pg.tiff", res = 600, units = "in", height = 5, width = 8)
sex_facet_hist
#dev.off()

# model validation ----
## validate trait type ----
#response = R2
#expl = Kingsolver_traits
#random = StudyID
# glm name = all.model.traits


# zuur book

#step 1 linear regression
lmTraits <- lm(R.2 ~ Kingsolver_traits, data = data.all.traits)
plot(lmTraits)
data.all.traits$logR.2 <- log10(data.all.traits$R.2 + 1)
loglmtraits <- lm(logR.2 ~ Kingsolver_traits, data = data.all.traits)
asinlmtraits <- lm(asin(R.2) ~ Kingsolver_traits, data = data.all.traits)
sqrtlmtraits <- lm(sqrt(R.2) ~ Kingsolver_traits, data = data.all.traits)


# visual inspection 
par(mfrow = c(2,2))
plot(sqrtlmtraits, add.smooth = FALSE, which = 1) # homogeneity (fitted values vs residuals)
hist(resid(sqrtlmtraits), xlab = "Residuals", main = "") # normality - not normal
plot(data.all.traits$Kingsolver_traits, resid(sqrtlmtraits), # note that spread not the same 
     xlab = "Trait type", ylab = "residuals")
par(op)

# non-visual test of homogeneity (bartlett)
# null hypothesis is that variances are equal 
bartlett.test(resid(lmTraits), data.all.traits$Kingsolver_traits) # reject null

#step 2 GLS
traitsForm <- formula(R.2 ~ Kingsolver_traits) 
glsTraits <- gls(traitsForm, data = data.all.traits)

#step 3 variance
# choose random effect?

#step 4 fit model
lmmTraits <- lme(traitsForm, random = ~ 1 | StudyID,
                 method = "REML", data = data.all.traits)

#step 5 compare new/old models
anova(glsTraits, lmmTraits) #model w random is better; L = 89.46 (df = 8, p < 0.0001)

#step6 is it ok?
residTraits <- resid(lmmTraits, type = "normalized")
fittedTraits <- fitted(lmmTraits)
par(mfrow = c(1,2), mar = c(4,4,3,2))
plot(x = fittedTraits, y = residTraits, xlab = "Fitted values", ylab = "Residuals")
boxplot(residTraits ~ Kingsolver_traits, data = data.all.traits,
        main = "Kingsolver traits", lab= "Residuals")
par(op)

#step7/8 optimal fixed str
# this does not apply because we have 1 fixed ??

#step9 validate the omodel
summary(lmmTraits)
# (0.188^2)/(0.188^2 + 0.202^2) = correlation between obs from the same stuyd = 0.46

#steo10

#step11
library(lattice)
xyplot(residTraits ~ Kingsolver_traits, 
       data = data.all.traits, 
       ylab = "Residuals",
       xlab = "Kingsolver_traits",
       panel = function(x,y)
         {panel.grid(h = -1, v = 2) 
         panel.points(x, y, col = 1) 
         panel.loess(x, y, span = 0.5, col = 1,lwd=2)
         }
       )


library(mgcv)

gammTraits <- gamm(R.2 ~ Kingsolver_traits,
                   random = list(StudyID = ~1),
                   data = data.all.traits) 
summary(gammTraits$gam)
summary(gammTraits$lme)

anova(gammTraits$gam)
anova(gammTraits$lme)

plot(gammTraits$gam, all.terms = TRUE)
plot(gammTraits$lme)

drop1(sex.and.triats)


# I have no idea what any of the above is showing atm 

# so, dharma package instead

library("DHARMa")

testDispersion(all.model.traits)
simulationOutput <- simulateResiduals(fittedModel = all.model.traits, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantilefunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)

testZeroInflation(simulationOutput)

countOnes <- function(x) sum(x == 1)
testGeneric(simulationOutput, summary = countOnes, alternative = "greater") 

M2 <- glmer.nb(R.2 ~ Kingsolver_traits + (1|StudyID),
             data = data.all.traits)

summary(M2)

E2 <- resid(M2, type = "pearson")
N  <- nrow(data.all.traits)
p  <- length(coef(M2)) + 1  # '+1' is for variance parameter in NB
sum(E2^2) / (N - p)


f = glmm.zinb(fixed = R.2 ~ Kingsolver_traits,
              random = ~ 1 | StudyID, data = data.all.traits) 
summary(f)
fixed(f)
summary(f$fit.zero)

## validate sex ----
all.model.sex

lmSex <- lm(R.2 ~ Sex, data = data.all)


par(mfrow = c(2,2))
plot(lmSex, add.smooth = FALSE, which = 1) # homogeneity (fitted values vs residuals)
hist(resid(lmSex), xlab = "Residuals", main = "") # normality - not normal
plot(data.all$Sex, resid(lmSex), # note that spread not the same 
     xlab = "Sex", ylab = "residuals")
par(op)

loglmSex <- lm(log10(R.2) ~ Sex, data = data.all)

plot(resid(all.model.sex))
hist(resid(all.model.sex))

## validate rearing ----
all.model.rearing
hist(resid(all.model.rearing))
plot(resid(all.model.rearing))
par(op)


# traits
plot(resid(all.model.traits))
hist(resid(all.model.traits))

# ecology
ecology.full

ecology.test <- glmer(R.2 ~ method + (1 | StudyID), 
                      family = binomial, data= data.for.ecology.models)
plot(resid(ecology.test))
hist(resid(ecology.test))
summary(ecology.full)

# evil
plot(resid(evolhist.full))
hist(resid(evolhist.full))
summary(evolhist.full)
# intro
# regression
data.for.intro.models.broad <- data.for.intro.models.broad %>% filter(!TraitID == 303)

introlm <- lm(R.2 ~ method, data = data.for.intro.models.broad)
plot(introlm) # trash
data.for.intro.models.broad$logR.2 <- log10(data.for.intro.models.broad$R.2 + 1)
introlmlog <- lm(logR.2 ~ method, data = data.for.intro.models.broad)
E <- rstandard(introlmlog)
plot(introlmlog)


# step2/3
introformula <- formula(logR.2 ~ method)
introgls <- gls(introformula, data = data.for.intro.models.broad)

# step4
introlme <- lme(introformula, random = ~ 1 | StudyID,
                method = "REML", data = data.for.intro.models.broad)

# step5
anova(introgls, introlme) # better with RE (L = 28.66, p < 0.0001)

#step6
normresidintrolme <- resid(introlme, type = "normalized")
fittedintrolme <- fitted(introlme)
plot(x = fittedintrolme, y = normresidintrolme)
boxplot(normresidintrolme ~ method,
        data = data.for.intro.models.broad)

#step7/8 - skip

#step9
summary(introlme)
# 0.66 = correlation of obs from same study

#step10 p 156
library(lattice)
xyplot(normresidintrolme ~ method,
                data = data.for.intro.models.broad,
       ylab = "residuals",
                panel = function(x,y){
                  panel.grid(h = -1, v = 2)
                  panel.points(x, y, col = 1)
                  panel.loess(x, y, span = 0.5, col = 1,lwd=2)})

library(mgcv)
introgamm <- gamm(logR.2 ~ method,
                  random = list(StudyID = ~ 1),
                  data = data.for.intro.models.broad)
summary(introgamm$gam)
summary(introgamm$lme)
anova(introgamm$gam) 
plot(introgamm$gam) # we didn't have smoother (only categorical data)
plot(introgamm$lme)

plot(resid(introlm))
hist(resid(introlm))
plot(resid(intro.full.broad))
hist(resid(intro.full.broad))


summary(intro.full.broad)




M5 <- glm(R.2 ~ method, family = binomial, data = data.for.intro.models.broad)

EP <- resid(M5, type = "pearson")
ED <- resid(M5, type = "deviance")
mu <- predict(M5, type = "response")
E <- data.for.intro.models.broad$R.2 ~ mu
par(mfrow = c(2,2))
plot(x = mu, y = E, main = "response residuals")
plot(x = mu, y = EP)
plot(x = mu, y = ED)
plot(M5)



test <- glm(R.2 ~ method, family = quasibinomial, data = data.for.intro.models.broad)
summary(test)
drop1(test, test = "F")


EP <- resid(test, type = "pearson")
ED <- resid(test, type = "deviance")
mu <- predict(test, type = "response")
E <- data.for.intro.models.broad$R.2 ~ mu
par(mfrow = c(2,2))
plot(x = mu, y = E, main = "response residuals")
plot(x = mu, y = EP)
plot(x = mu, y = ED)
plot(test)



a1b <- (glmer.nb(R.2 ~ method + (1 | StudyID), data = data.for.intro.models.broad))
b2b <- (glmer.nb(R.2 ~ method + (1 | StudyID), data = data.for.ecology.models)) 
c2b <- (glmer.nb(R.2 ~ method + (1 | StudyID), data = data.for.evolhist.models)) 



library("DHARMa")

testDispersion(c2b)
simulationOutput <- simulateResiduals(fittedModel = c2b, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantilefunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)

testZeroInflation(simulationOutput)

countOnes <- function(x) sum(x == 1)
testGeneric(simulationOutput, summary = countOnes, alternative = "greater") 

library(NBZIMM)
f = glmm.zinb(fixed = log10(R.2+1) ~ method, 
              random = ~ 1 | StudyID, data = data.for.intro.models.broad) 
summary(f)
plot(f)
plot(intro.full.broad)

g = glmm.zinb(fixed = R.2 ~ method, 
              random = ~ 1 | StudyID, data = data.for.ecology.models) 
summary(g)
plot(g)
plot(ecology.full)

h = glmm.zinb(fixed = R.2 ~ method, 
                random = ~ 1 | StudyID, data = data.for.evolhist.models) 
summary(h)
plot(h)
plot(evolhist.full)

plot(resid(h))
plot(resid(f))
plot(resid(h))
plot(resid(introlm))
hist(resid(introlm))
hist(resid(h))
hist(resid(g))
hist(resid(f))
plot(resid(intro.full.broad))
hist(resid(intro.full.broad))

all.model.rearing

testDispersion(all.model.rearing)
simulationOutput <- simulateResiduals(fittedModel = all.model.rearing, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantilefunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)

testZeroInflation(simulationOutput)


## ecology

ecology.full

testDispersion(ecology.full)
simulationOutput <- simulateResiduals(fittedModel = ecology.full, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantilefunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)

testZeroInflation(simulationOutput)



## evilhist

evolhist.full

testDispersion(evolhist.full)
simulationOutput <- simulateResiduals(fittedModel = evolhist.full, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantilefunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)

testZeroInflation(simulationOutput)



## introductions

intro.full.broad

testDispersion(intro.full.broad)
simulationOutput <- simulateResiduals(fittedModel = intro.full.broad, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantilefunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
plotQQunif(simulationOutput)
plotResiduals(simulationOutput)

testZeroInflation(simulationOutput)

