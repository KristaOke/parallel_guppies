library(stats)
library(tidyverse)

#setwd("<YOURPATHWAY>")
wd <- getwd()
new.data<- read.csv(paste(wd,'/Data/MetaData.csv',sep=""), header=TRUE, sep=",")

#my data was read in as a csv file with columns: 

# "Study.ID"             "Species"              "Contrast"             "Trait.ID"             "Trait"                "Trait.type"           "Standardized."       
# "Standardization.type" "Paired."              "Study.type"           "Data.point"           "Treatment"            "Population"           "Mean"                
# "N" 

#select only entries with both predation levels
a<-new.data
a$High<- ifelse(new.data$Predation== "High", print(new.data$TraitID), NA )
a$Low<- ifelse(new.data$Predation== "Low", print(new.data$TraitID), NA )

a$High[!(a$High %in% a$Low)] #none feb 10th
a$Low[!(a$Low %in% a$High)] #none

test<-new.data %>% 
  group_by(TraitID, Predation) %>% 
  tally() %>% 
  filter(n<2) #looks for Traits with less than 4 (2 high and 2 low) populations, none feb 10th


#filter out incomplete entries, and entries with only one data point for either high or low (will make R2 1.0)
new.data<-new.data  %>% 
  filter(!is.na(Predation)) %>% 
  filter(!is.na(MeanValue)) 
  #filter(!(TraitID %in% c(??))) #these are entries that need to be excluded (only one pred level), none feb 10


new.data$Predation<-as.factor(new.data$Predation)
new.data$MeanValue<-as.numeric(new.data$MeanValue)

#how many traits/studies?
length(unique(new.data$TraitID))#555 for TraitID (without typos traits) jul 21
length(unique(new.data$StudyID)) #48

#set up a new function to loop through all data and run anova on Mean trait values for each trait
#note this assumes we have just one covariate of interest: predation(treatment) (basically habitat, e.g. high/low predation)
#if we want more covariates we'll need to rethink this a bit

#START FUNC

dat<-FALSE 

i <- 1

anova_loop <- function(dat=dat){
  
  trait <- as.factor(dat$TraitID)
  
  output <- matrix(ncol=7,nrow=length(levels(trait)))
  colnames(output) <- c("R^2","adj.R^2","partial.Eta","F","p","TraitID", "justr")
  
  for (i in 1: length(levels(trait))){
    
    print(i)
    
    aov.mod   <- aov(MeanValue~Predation,data=dat[dat$TraitID==levels(trait)[i],])
    anova.mod <- stats:::anova.lm(aov.mod)
    mod.ss    <- anova.mod$"Sum Sq"
    mod.pes   <- mod.ss/(mod.ss+mod.ss[length(mod.ss)])          # calculate pes
    
    r <- summary.lm(aov.mod)                                     # summary of linear model
    
    R     <- round(r$"r.squared",3)                              # R squared
    R_adj <- round(r$"adj.r.squared",3)                          # Adjusted R squared
    F_value <- round(as.numeric(anova.mod$"F value")[1],3)       # F
    p_value <- round(as.numeric(anova.mod$"Pr(>F)")[1],3)        # p value
    partial_eta <- round(mod.pes[1],3)                           # pes
    TraitID <- levels(trait)[i]
    justr <- sqrt(anova.mod[1,4]/(anova.mod[1,4]+anova.mod[2,1]))
    
    output[i,] <- c(R,R_adj,partial_eta,F_value,p_value, TraitID, justr) 
  }
  output <- data.frame(output)
  return(output)
  
}

#END FUNC


output.all<-anova_loop(dat=new.data) #use function to run loop on all data

#look at some histograms of output
hist(as.numeric(as.character(output.all$partial.Eta)))
hist(as.numeric(as.character(output.all$R.2)))
hist(as.numeric(as.character(output.all$adj.R.2)))
hist(as.numeric(as.character(output.all$F)))
hist(as.numeric(as.character(output.all$justr)))


#write results to a csv file for further analyses DONT FORGET TO MOVE AND OVERWRITE saves inside repo not data folder

#this is the overall R2, all studies and traits included
#write.table(output.all, file = "TraitR2_among.csv",row.names=FALSE,col.names=TRUE, sep=",") 

# anova for the south only (Slope Q)
output.south<-anova_loop(dat=traits_s)
#write.table(output.south, file = "TraitR2_South.csv",row.names=FALSE,col.names=TRUE, sep=",") 

# anova for the caroni only (drainage Q)
output.caroni<-anova_loop(dat=traits_d)
#write.table(output.caroni, file = "TraitR2_Caroni.csv",row.names=FALSE,col.names=TRUE, sep=",") 

# anova for the south only (drainage Q, northern pops excluded because they arent in caroni or oropuche drainages)
output.south.drainages<-anova_loop(dat=traits_d_all)
#write.table(output.south.drainages, file = "TraitR2_Among_Drainage.csv",row.names=FALSE,col.names=TRUE, sep=",") 

# anova for the natural pops only (drainage Q)
output.intros<-anova_loop(dat=traits_i)
#write.table(output.intros, file = "TraitR2_intro.csv",row.names=FALSE,col.names=TRUE, sep=",")

# anova for the natural pops only (drainage Q, broad specification)
output.intros.broad<-anova_loop(dat=traits_i_broad)
#write.table(output.intros, file = "TraitR2_intro_broad.csv",row.names=FALSE,col.names=TRUE, sep=",") 