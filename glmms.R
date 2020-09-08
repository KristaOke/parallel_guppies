# GLMMs for parallel_guppies! 
# 2020-08-31 AH
# ok

# I've kept in all of the models that don't work, but you can DELETE THEM.
# Just wanted them in so you could see what I've tried so far 

# LIBRARIES ---- 
library(dplyr)
library(tidyr)
library(lme4)
library(ggplot2)
library(bbmle)
library(gridExtra)
library(car)
library(corrplot)
library(tidyverse)

# set working directory
setwd("~/R/MSc/Side_projects/Meta_guppies_parallel/Data")

# Import and tidy
spreadsheet.data <- read.csv("MetaData_sexID.csv")
R2.data <- read.csv("TraitR2_sexID.csv")

# spreadsheet.data is the data extracted for the meta-analysis
# R2.data is the output of the ANOVA loop

names(R2.data)[names(R2.data) == "TraitID"] <- "sex_TraitID"  # so same in both spreadsheets

str(spreadsheet.data)

spreadsheet.data$Study.ID <- as.factor(spreadsheet.data$Study.ID)
spreadsheet.data$Collection_start <- as.factor(spreadsheet.data$Collection_start)
spreadsheet.data$Collection_end <- as.factor(spreadsheet.data$Collection_end)
spreadsheet.data$Published <- as.factor(spreadsheet.data$Published)
spreadsheet.data$TraitID <- as.factor(spreadsheet.data$TraitID)

str(R2.data)
R2.data$sex_TraitID <- as.factor(R2.data$sex_TraitID)

# This (data.for.models) is the data to use
data.for.models <- merge(spreadsheet.data, R2.data,  by = "sex_TraitID")
data.for.models$Sex <- gsub("Both (mostly juveniles)",  # this is just for now 
                            "Both", 
                            data.for.models$Sex, 
                            fixed = T)  # fixed in spreadsheet

data.for.models$Sex <- as.factor(data.for.models$Sex)
data.for.models$sex_TraitID <- as.factor(data.for.models$sex_TraitID)

str(data.for.models)

##%######################################################%##
#                                                          #
####                  EXPLORING DATA                    ####
#                                                          #
##%######################################################%##


plot(table(data.for.models$R.2), ylab = "Frequency", xlab = "R.2")

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
  group_by(Slope, Sex, TraitType2, Paired) %>% tally())

(wctally <-
  data.for.models %>% 
  filter(StudyType == "Wildcaught" & Slope %in% c("North", "South") & Sex %in% c("M", "F", "Both")) %>% 
  group_by(Slope, Sex, TraitType2, Paired) %>% tally())
  
## MOD 1 Traits/both sexes

# WC 
(overall.trait.plot.WC <-  # many more studies in the south....  
    data.for.models %>% filter(StudyType == "Wildcaught" & Sex == "Both") %>% 
    ggplot(aes(x = TraitType2, y = R.2)) +
    theme_bw() + 
    geom_violin(aes(fill = TraitType2)) + 
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

tbl1 <- table(data.for.models$TraitType2, data.for.models$Slope)
chisq.test(tbl1)  

tbl2 <- table(data.for.models$TraitType2, data.for.models$Paired)
chisq.test(tbl2)  

tbl3 <- table(data.for.models$Paired, data.for.models$Slope)
chisq.test(tbl3)  

corrplot(tbl1, is.cor = FALSE)
corrplot(tbl2, is.cor = FALSE)
corrplot(tbl3, is.cor = FALSE)

tbl4 <- table(data.for.models$TraitType2, data.for.models$Sex)
tbl4
corrplot(tbl4, is.cor = FALSE)


##%######################################################%##
#                                                          #
####                  MODELS                            ####
#                                                          #
##%######################################################%##


#### QUESTION 1. IS THERE A DIFFERENCE OVERALL BETWEEN HIGH/LOW ####
## Traittype2 is the fixed effect that we decided  we really care about here ## 

# I think this one works THE BEST -- BUT only South
mod1.south <- glmer(R.2 ~ TraitType2 + (1|Study.ID), 
                       family = binomial,
              data = data.for.models[
                data.for.models$StudyType == "Wildcaught" 
                & data.for.models$Sex == "Both"
                & data.for.models$Slope == "South",
                ])
summary(mod1.south)
AIC(mod1.south)  # 494.1767

# has drainages, won't converge - crap
mod1.drain <- glmer(R.2 ~ TraitType2 + (1|Study.ID/Drainage), 
                    family = binomial, 
                    data = data.for.models[
                      data.for.models$StudyType == "Wildcaught" 
                      & data.for.models$Sex == "Both",])

summary(mod1.drain)
AIC(mod1.drain)

# Again, this one won't converge - crappy crap
mod1.slope <- glmer(R.2 ~ TraitType2 + (1|Study.ID/Slope), 
                    family = binomial, 
                    data = data.for.models[
                      data.for.models$StudyType == "Wildcaught" 
                      & data.for.models$Sex == "Both",])

summary(mod1.slope)
AIC(mod1.slope)

# Only North - Won't run (error) 
mod1.north <- glmer(R.2 ~ TraitType2 + (1|Study.ID), 
                           family = binomial,
                           data = data.for.models[
                             data.for.models$StudyType == "Wildcaught" 
                             & data.for.models$Sex == "Both"
                             & data.for.models$Slope == "North",
                             ])

# This follows totally to our original plan... again, model won't converge
trait.mod1 <- glmer(R.2 ~ Slope + TraitType2 + Paired + (1|Study.ID), 
                           family = binomial,
                           data = data.for.models[
                             data.for.models$StudyType == "Wildcaught" 
                             & data.for.models$Sex == "Both",
                             ])

summary(trait.mod1)
AIC(trait.mod1)  # 551.0687

# Runs but singular fit... TraitType 2 is mostly Paired anyway... crap
trait.mod2 <- glmer(R.2 ~ TraitType2 + Paired + (1|Study.ID), 
                       family = binomial,
                       data = data.for.models[
                         data.for.models$StudyType == "Wildcaught"
                         & data.for.models$Sex == "Both",
                         ])

summary(trait.mod2)
AIC(trait.mod2)  # 591.8565

#### QUESTION 2. IS PARALLELISM DIFFERENT BETWEEN THE SEXES? ####

sex.mod1 <- glmer(R.2 ~ Sex*Slope + (1|Study.ID),
                   family = binomial, 
                   data = data.for.models[data.for.models$StudyType == "Wildcaught"
                                          & data.for.models$Sex %in% c("M", "F"),])
summary(sex.mod1)
AIC(sex.mod1) # 1995.517

sex.mod2 <- glmer(R.2 ~ Sex+Slope + (1|Study.ID),
                     family = binomial, 
                     data = data.for.models[data.for.models$StudyType == "Wildcaught"
                                            & data.for.models$Sex %in% c("M", "F"),])
summary(sex.mod2)
AIC(sex.mod2)  # 1994.825

sex.mod3 <- glmer(R.2 ~ Sex + (1|Study.ID),
                     family = binomial, 
                     data = data.for.models[data.for.models$StudyType == "Wildcaught"
                                            & data.for.models$Sex %in% c("M", "F"),])
summary(sex.mod3)
AIC(sex.mod3)  # 2387.036


# This dude won't coverge
sex.mod4 <- glmer(R.2 ~ Sex*TraitType2 + (1|Study.ID),
                     family = binomial, 
                     data =
                       data.for.models[data.for.models$StudyType == "Wildcaught"
                                       & data.for.models$Sex %in% c("M", "F")
                                       & data.for.models$TraitType2 %in% c("Diet", "Morphometric", "Other", "Physiology", "Behaviour"),])
summary(sex.mod4)
AIC(sex.mod4)

sex.mod5 <- glmer(R.2 ~ Sex + TraitType2 + (1|Study.ID),
                  family = binomial, 
                  data =
                    data.for.models[data.for.models$StudyType == "Wildcaught"
                                    & data.for.models$Sex %in% c("M", "F")
                                    & data.for.models$TraitType2 %in% c("Diet", "Morphometric", "Other", "Physiology", "Behaviour"),])
summary(sex.mod5)
AIC(sex.mod5)

###


sex.mod7.wc <- glmer(R.2 ~ Sex + (1|Study.ID/Slope),
                     family = binomial, 
                     data =
                       data.for.models[data.for.models$StudyType == "Wildcaught"
                                       & data.for.models$Sex %in% c("M", "F")
                                       & data.for.models$TraitType2 %in% c("Diet", "Morphometric", "Other", "Physiology", "Behaviour"),])
summary(sex.mod7.wc)
AIC(sex.mod7.wc)  # 1190/876

# Males

trait.mod1.males <- glmer(R.2 ~ Slope + TraitType2 + Paired + (1|Study.ID), 
                          family = binomial,
                  data = data.for.models[
                    data.for.models$StudyType == "Wildcaught" 
                    & data.for.models$Sex == "M",
                    ])

summary(trait.mod1.males)
AIC(trait.mod1.males)  # -1972.583

trait.mod2.males <- glmer(R.2 ~ Slope + Paired + (1|Study.ID), 
                          family = binomial,
                          data = data.for.models[data.for.models$StudyType == "Wildcaught" & data.for.models$Sex == "M",])

summary(trait.mod2.males)
AIC(trait.mod2.males)  # 1307.4

trait.mod3.males <- glmer(R.2 ~ Slope + TraitType2 + (1|Study.ID), 
                          family = binomial,
                          data = data.for.models[
                            data.for.models$StudyType == "Wildcaught" 
                            & data.for.models$Sex == "M",
                            ])

summary(trait.mod3.males)
AIC(trait.mod3.males)  # 1292.874

trait.mod4.males <- glmer(R.2 ~ TraitType2 + Paired + (1|Study.ID), 
                          family = binomial,
                          data = data.for.models[
                            data.for.models$StudyType == "Wildcaught" 
                            & data.for.models$Sex == "M",
                            ])

summary(trait.mod4.males)
AIC(trait.mod4.males)  # 1497.781

trait.mod5.males <- glmer(R.2 ~ TraitType2 + (1|Study.ID), 
                          family = binomial,
                          data = data.for.models[
                            data.for.models$StudyType == "Wildcaught" 
                            & data.for.models$Sex == "M",
                            ])

summary(trait.mod5.males)
AIC(trait.mod5.males)  # 1495.068

trait.mod6.males <- glmer(R.2 ~ TraitType2 + (1|Study.ID), 
                          family = binomial,
                          data = data.for.models[
                            data.for.models$StudyType == "Wildcaught" 
                            & data.for.models$Sex == "M"
                            & data.for.models$Slope == "South",
                            ])

summary(trait.mod6.males)
AIC(trait.mod6.males)  # 806.9988

trait.mod7.males <- glmer(R.2 ~ TraitType2 + (1|Study.ID), 
                          family = binomial,
                          data = data.for.models[
                            data.for.models$StudyType == "Wildcaught" 
                            & data.for.models$Sex == "M"
                            & data.for.models$Slope == "North",
                            ])

summary(trait.mod7.males)
AIC(trait.mod7.males)  # 806.9988



# Females

trait.mod1.females <- glmer(R.2 ~ Slope + TraitType2 + Paired + (1|Study.ID), 
                            family = binomial,
                            data = data.for.models[
                              data.for.models$StudyType == "Wildcaught" 
                              & data.for.models$Sex == "F",
                              ])

summary(trait.mod1.females)
AIC(trait.mod1.females)  # 636.5785

#### QUESTION 3 - IS THERE A DIFFERENCE BETWEEN THE SLOPES? ####

slope.mod1 <- glmer(R.2 ~ Slope + (1|Study.ID),
                     family = binomial, 
                     data = data.for.models[data.for.models$StudyType == "Wildcaught"
                                            & data.for.models$Sex == "Both",])
summary(slope.mod1)
AIC(slope.mod1)  # 599.0427

#### QUESTION 4. IS THERE REGRESSION TOWARDS THE MEAN ####
## Here, we are looking at both sexes and Collection_end ##

# Honestly Allegra I don't know what to do about this mess below...
# Will also have to separate by TraitType too, which I forgot about here.

data.for.models %>% 
  filter(Collection_end != "NA") %>% 
  group_by(Study.ID, Collection_end) %>% 
  ggplot(aes(x = R.2)) + geom_histogram(binwidth = 0.03) 

data.reg2m <- data.for.models %>% 
  filter(Collection_end !="NA") %>% 
  group_by(Collection_end, Study.ID) %>% 
  summarise(meanR2 = mean(R.2))

head(data.reg2m)

# individual R2 for Collection_end years 
# don't know what to do because not paired ... 

R2 <- data.reg2m$meanR2
names(R2) <- paste(data.reg2m$Collection_end)
print(R2)
mean(R2)  

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
  geom_hline(yintercept = 0.3770349, colour = "blue", linetype = "dashed")

# These notes below are how you would calculate regression to the meet if you had two time points per id
# From Mee and Chua - regression towards the mean and the paired sample t-test (1991)  

# 1. Calculate X = Y1-??
# 2. Estimate the parameters ?? 0 and ?? from the linear regression model of Y2 on X
# 3. Estimate the treatment effect ??^ by subtracting ?? from ??^0, the estimate of ?? 0
# 4. Calculate the test-statistic
