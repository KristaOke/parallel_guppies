# Analysis for 'Compiling forty years of guppy research to investigate the factors contributing to (non)parallel evolution'
#This R script uses the prepared R2 files provided on the repository, see "GuppyR2_Prep.R" for the process of generating these

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


# import and tidy R2 values to prepare for running models (descriptions of documents in ReadMe)

spreadsheet.data <- read.csv(paste(wd,'/Data/MetaData.csv',sep=""), header=TRUE, sep=",")
R2.data.among <- read.csv(paste(wd,'/Data/TraitR2_among.csv',sep=""), header=TRUE, sep=",")
R2.data.south <- read.csv(paste(wd,'/Data/TraitR2_south.csv',sep=""), header=TRUE, sep=",")
#R2.data.intro <- read.csv(paste(wd,'/Data/TraitR2_intro.csv',sep=""), header=TRUE, sep=",") # Not used in paper (doesn't include all introduction information in Table S2)
R2.data.caroni <- read.csv(paste(wd,'/Data/TraitR2_Caroni.csv',sep=""), header=TRUE, sep=",")
R2.data.among.drainage <- read.csv(paste(wd,'/Data/TraitR2_Among_Drainage.csv',sep=""), header=TRUE, sep=",")
R2.data.intro.broad <- read.csv(paste(wd, '/Data/TraitR2_intro_broad.csv', sep = ""), header = TRUE, sep = ",")

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
R2.data.intro.broad$TraitID <- as.factor(R2.data.intro.broad$TraitID)
R2.data.caroni$TraitID <- as.factor(R2.data.caroni$TraitID)
R2.data.among.drainage$TraitID <- as.factor(R2.data.among.drainage$TraitID)

# combine R2 and spreadsheet data
## prep for when we bind them together, this variable will go in the model to indicate the type of R2
R2.data.among$method <- "all"
R2.data.south$method <- "south"
R2.data.intro.broad$method <- "only_natural_broad"
R2.data.caroni$method <- "caroni"
R2.data.among.drainage$method <- "both.drainages"

## get relevant data with one entry for each Trait
data.all <- inner_join(spreadsheet.data, R2.data.among, by = "TraitID") %>% 
  dplyr::select(1:3, 6:14, 17:21, 23, 43:52)%>% #selected only the columns that have information that applies at the trait level
  distinct(TraitID, .keep_all = TRUE) #removes duplicated TraitID rows and retain columns

## repeat for other subset R2
data.south <- inner_join(spreadsheet.data, R2.data.south, by = "TraitID") %>% 
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

### Models for basic potential determinants of parallelism -----
# remove 'both' from sex (only used traits measured in male or female guppies)
(data.all <- data.all %>% filter(Sex %in% c("M", "F"))) %>% mutate(Sex = droplevels(Sex))

## Trait type -------
## remove 'other'
data.all.no.other <- data.all %>% filter(!Kingsolver_traits == "Other") %>% mutate(Kingsolver_traits = droplevels(Kingsolver_traits)) 

(all.model.traits <- glmer(R.2 ~ Kingsolver_traits +  (1|StudyID), 
                           data = data.all.no.other, family = binomial)) %>% summary()

car::Anova(all.model.traits, type = "II") # anova to get Chi-sq
## Sex (with and without colour) -------
# sex with colour
(all.model.sex <- glmer(R.2 ~ Sex + (1|StudyID), data = data.all, family = binomial)) %>% summary()

car::Anova(all.model.sex, type = "II") # anova to get Chi-sq

# sex without colour
# remove 'colour'
(data.all.no.colour <- data.all %>% filter(!Kingsolver_traits == 'Colour')) %>% mutate(Kingsolver_traits = droplevels(Kingsolver_traits))

(all.model.sex.no.colour <- glmer(R.2 ~ Sex + (1|StudyID), data = data.all.no.colour, family = binomial)) %>% summary()

car::Anova(all.model.sex.no.colour, type = "II") # anova to get Chi-sq
## Rearing environment ------
data.all.rear <- data.all %>% filter(StudyType %in% c("Common Garden (F2)", "Wildcaught"))

(all.model.rearing <- glmer(R.2 ~ StudyType +  (1|StudyID), data = data.all.rear, family = binomial)) %>% summary()

car::Anova(all.model.rearing, type = "II") # anova to get Chi-sq

## Sex and rearing environment ------
# Use the same df as in rearing environment
(sex.and.rear <- glmer(R.2 ~ StudyType + Sex + (1|StudyID), data = data.all.rear, family = binomial)) %>% summary()

Anova(sex.and.rear, type = 2) # anova to get Chi-sq

### Ecological Complexity (North vs Both slopes) -----
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

### Evolutionary History (Caroni vs entire South slope) ----
#Note: Northern slope excluded from this subset
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

### Time since Colonization (Natural vs Introduced and natural) ----
data.for.intro.models.broad<-rbind(data.all,data.intro.broad) %>% 
  arrange(TraitID) %>% #puts them in a nice order
  group_by(TraitID) %>% #groups them for the count
  filter(n() > 1) #filters only trait IDs that have more than 1 entry within the group (n = 204 traits)

### Permutation test ----





### Other values (Mean and SD) reported in text ----

count(data.all %>% filter(R.2 < 0.5))  # 364
count(data.all)  # 446
364/446  # proportion of R2 values less than 0.5

count(data.all %>% filter(R.2 < 0.1))
190/446

count(data.all %>% filter(R.2 > 0.9))
6/446
# Mean R2
mean(data.all$R.2)
sd(data.all$R.2)

## Traits
(Colourtrait <- data.all %>% filter(Kingsolver_traits == "Colour")) %>% summary()
sd(Colourtrait$R.2)
(LHtrait <- data.all %>% filter(Kingsolver_traits == "Other_life_history")) %>% summary()
sd(LHtrait$R.2)
(Sizetrait <- data.all %>% filter(Kingsolver_traits == "Size")) %>% summary()
sd(Sizetrait$R.2)
(Morphtrait <- data.all %>% filter(Kingsolver_traits == "Other_morphology")) %>% summary()
sd(Morphtrait$R.2)
(Othertrait <- data.all %>% filter(Kingsolver_traits == "Other")) %>% summary()
sd(Othertrait$R.2)
(Phystrait <- data.all %>% filter(Kingsolver_traits == "Physiology")) %>% summary()
sd(Phystrait$R.2)
(Behavtrait <- data.all %>% filter(Kingsolver_traits == "Behaviour")) %>% summary()
sd(Behavtrait$R.2)

## Sex
(Malewithcolour <- data.all %>% filter(Sex == "M")) %>% summary()
sd(Malewithcolour$R.2)
(Malenocolour <- data.all.no.colour %>% filter(Sex == "M")) %>% summary()
sd(Malenocolour$R.2)
(Femalesex <- data.all %>% filter(Sex == "F")) %>% summary()
sd(Femalesex$R.2)

## Rearing enviro
(cg <- data.all %>% filter(StudyType == "Common Garden (F2)")) %>% summary()
sd(cg$R.2)
(wc <- data.all %>% filter(StudyType == "Wildcaught")) %>% summary()
sd(wc$R.2)

## Ecology
(one.slope <- data.for.ecology.models %>% filter(method == "south")) %>% summary()
sd(one.slope$R.2)
(both.slopes <- data.for.ecology.models %>% filter(method == "all")) %>% summary()
sd(both.slopes$R.2)

## Intro
data.for.intro.models.broad <- data.for.intro.models.broad %>% drop_na(R.2)
(only_natural_broad <- data.for.intro.models.broad %>% filter(method == "only_natural_broad")) %>% summary()
sd(only_natural_broad$R.2)
(natural_and_intro_broad <- data.for.intro.models.broad %>% filter(method == "all")) %>% summary()
sd(natural_and_intro_broad$R.2)

## Evolhist
(one_drainage <- data.for.evolhist.models %>% filter(method == "caroni")) %>% summary()
sd(one_drainage$R.2)

(both_drainages <- data.for.evolhist.models %>% filter(method == "both.drainages")) %>% summary()
sd(both_drainages$R.2)


### figures in manuscript ----

## figure 1 is map 

## figure 2
fig2data <- read.csv("testPG_figure.csv", fileEncoding="UTF-8-BOM")
fig2data$population <- as.factor(fig2data$population)

fig2data$population <- factor(fig2data$population, levels = (c("Marianne high", 
                                                               "Marianne low", 
                                                               "Aripo high", "Aripo low",
                                                               "Damier high",
                                                               "Damier low",
                                                               "Yarra high",
                                                               "Yarra low")))

levels(fig2data$population)

highparal <- fig2data %>% filter(facet == 'a') %>% 
  ggplot(aes(x = population, y = mean, color = predation)) +
  geom_point(size = 4) +
  labs(title = "Example of high R  (0.984)",
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
        plot.title = element_text(size = 16),
        axis.text.x = element_text(angle = 45,  hjust = 1))

lowparal <- fig2data %>% filter(facet == 'b') %>% 
  ggplot(aes(x = population, y = mean, color = predation)) +
  geom_point(size = 4) +
  labs(title = "Example of low R  (0.000)",
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
        plot.title = element_text(size = 16),
        axis.text.x = element_text(angle = 45,  hjust = 1))

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

## figure 3 

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

## figure 4 
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

## figure 5 

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
# tiff("pg.fig5.tiff", res= 600, units = "in", height = 6, width =8)
fig5
dev.off()

## figure 6 

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


### Supplementary material ----
## Supplemental table 1 

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

## Supplemental table 3 

R2Female <- data.all %>% filter(Sex == "F")
mean(R2Female$R.2)
sd(R2Female$R.2)

mean(R2Male.noColour$R.2)
sd(R2Male.noColour$R.2)

R2Female.Colour <- data.all %>% filter(Sex == "F" &
                                         Kingsolver_traits == "Colour")
mean(R2Female.Colour$R.2)
sd(R2Female.Colour$R.2)

R2Female.Behaviour <- data.all %>% filter(Sex == "F" &
                                            Kingsolver_traits == "Behaviour")
mean(R2Female.Behaviour$R.2)
sd(R2Female.Behaviour$R.2)

R2Female.LH <- data.all %>% filter(Sex == "F" &
                                     Kingsolver_traits == "Other_life_history")
mean(R2Female.LH$R.2)
sd(R2Female.LH$R.2)

R2Female.Size <- data.all %>% filter(Sex == "F" &
                                       Kingsolver_traits == "Size")
mean(R2Female.Size$R.2)
sd(R2Female.Size$R.2)

R2Female.Morph <- data.all %>% filter(Sex == "F" &
                                        Kingsolver_traits == "Other_morphology")
mean(R2Female.Morph$R.2)
sd(R2Female.Morph$R.2)

R2Female.Physiology <- data.all %>% filter(Sex == "F" &
                                             Kingsolver_traits == "Physiology")
mean(R2Female.Physiology$R.2)
sd(R2Female.Physiology$R.2)

R2Female.Other <- data.all %>% filter(Sex == "F" &
                                        Kingsolver_traits == "Other")
mean(R2Female.Other$R.2)
sd(R2Female.Other$R.2)

R2Female.cg <- data.all.rear %>% filter(Sex == "F" &
                                          StudyType == c("Common Garden (F2)"))
mean(R2Female.cg$R.2)
sd(R2Female.cg$R.2)

R2Female.wc <- data.all.rear %>% filter(Sex == "F" &
                                          StudyType == c("Wildcaught"))
mean(R2Female.wc$R.2)
sd(R2Female.wc$R.2)


R2femaleeco.within <- data.for.ecology.models %>% filter(Sex == "F" &
                                                           method == "south")

mean(R2femaleeco.within$R.2)
sd(R2femaleeco.within$R.2)

R2femaleeco.between <- data.for.ecology.models %>% filter(Sex == "F" &
                                                            method == "all")

mean(R2femaleeco.between$R.2)
sd(R2femaleeco.between$R.2)

R2femaleevo.within <- data.for.evolhist.models %>% filter(Sex == "F" &
                                                            method == "caroni")

mean(R2femaleevo.within$R.2)
sd(R2femaleevo.within$R.2)

R2femaleevo.between <- data.for.evolhist.models %>% filter(Sex == "F" &
                                                             method == "both.drainages")

mean(R2femaleevo.between$R.2)
sd(R2femaleevo.between$R.2)

R2femaletime.within <- data.for.intro.models.broad %>% filter(Sex == "F" &
                                                                method == "only_natural_broad")

mean(R2femaletime.within$R.2)
sd(R2femaletime.within$R.2)

R2femaletime.between <- data.for.intro.models.broad %>% filter(Sex == "F" &
                                                                 method == "all")

mean(R2femaletime.between$R.2)
sd(R2femaletime.between$R.2)


R2Male.wcolour <- data.all %>% filter(Sex == "M")
mean(R2Male.wcolour$R.2)
sd(R2Male.wcolour$R.2)

R2Male.noColour <- data.all %>% filter(Sex == "M" &
                                         !Kingsolver_traits == "Colour")
mean(R2Male.noColour$R.2)
sd(R2Male.noColour$R.2)

R2Male.Colour <- data.all %>% filter(Sex == "M" &
                                       Kingsolver_traits == "Colour")
mean(R2Male.Colour$R.2)
sd(R2Male.Colour$R.2)

R2Male.Behaviour <- data.all %>% filter(Sex == "M" &
                                          Kingsolver_traits == "Behaviour")
mean(R2Male.Behaviour$R.2)
sd(R2Male.Behaviour$R.2)

R2Male.LH <- data.all %>% filter(Sex == "M" &
                                   Kingsolver_traits == "Other_life_history")
mean(R2Male.LH$R.2)
sd(R2Male.LH$R.2)

R2Male.Size <- data.all %>% filter(Sex == "M" &
                                     Kingsolver_traits == "Size")
mean(R2Male.Size$R.2)
sd(R2Male.Size$R.2)

R2Male.Morph <- data.all %>% filter(Sex == "M" &
                                      Kingsolver_traits == "Other_morphology")
mean(R2Male.Morph$R.2)
sd(R2Male.Morph$R.2)

R2Male.Physiology <- data.all %>% filter(Sex == "M" &
                                           Kingsolver_traits == "Physiology")
mean(R2Male.Physiology$R.2)
sd(R2Male.Physiology$R.2)

R2Male.Other <- data.all %>% filter(Sex == "M" &
                                      Kingsolver_traits == "Other")
mean(R2Male.Other$R.2)
sd(R2Male.Other$R.2)

R2Male.cg <- data.all.rear %>% filter(Sex == "M" &
                                        StudyType == "Common Garden (F2)")
mean(R2Male.cg$R.2)
sd(R2Male.cg$R.2)

R2Male.wc <- data.all.rear %>% filter(Sex == "M" &
                                        StudyType == "Wildcaught")
mean(R2Male.wc$R.2)
sd(R2Male.wc$R.2)

R2Maleeco.within <- data.for.ecology.models %>% filter(Sex == "M" &
                                                         method == "south")

mean(R2Maleeco.within$R.2)
sd(R2Maleeco.within$R.2)

R2Maleeco.between <- data.for.ecology.models %>% filter(Sex == "M" &
                                                          method == "all")

mean(R2Maleeco.between$R.2)
sd(R2Maleeco.between$R.2)

R2Maleevo.within <- data.for.evolhist.models %>% filter(Sex == "M" &
                                                          method == "caroni")

mean(R2Maleevo.within$R.2)
sd(R2Maleevo.within$R.2)

R2Maleevo.between <- data.for.evolhist.models %>% filter(Sex == "M" &
                                                           method == "both.drainages")

mean(R2Maleevo.between$R.2)
sd(R2Maleevo.between$R.2)

R2Maletime.within <- data.for.intro.models.broad %>% filter(Sex == "M" &
                                                              method == "only_natural_broad")

mean(R2Maletime.within$R.2)
sd(R2Maletime.within$R.2)

R2Maletime.between <- data.for.intro.models.broad %>% filter(Sex == "M" &
                                                               method == "all")

mean(R2Maletime.between$R.2)
sd(R2Maletime.between$R.2)


## Supplemental figure 1 

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
###