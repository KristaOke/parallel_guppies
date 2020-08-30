library(stats)

#setwd("<YOURPATHWAY>")
wd <- getwd()
new.data<- read.csv(paste(wd,'/data/MetaData.csv',sep=""), header=TRUE, sep=",")

#my data was read in as a csv file with columns: 

# "Study.ID"             "Species"              "Contrast"             "Trait.ID"             "Trait"                "Trait.type"           "Standardized."       
# "Standardization.type" "Paired."              "Study.type"           "Data.point"           "Treatment"            "Population"           "Mean"                
# "N" 


#select only entries with both predation levels
a<-new.data
a$High<- ifelse(new.data$Treatment== "High", print(new.data$TraitID), NA )
a$Low<- ifelse(new.data$Treatment== "Low", print(new.data$TraitID), NA )

a$High[!(a$High %in% a$Low)] #4
a$Low[!(a$Low %in% a$High)] #630

#sex_TraitID column info
a<-new.data
a$High<- ifelse(new.data$Treatment== "High", print(new.data$sex_TraitID), NA )
a$Low<- ifelse(new.data$Treatment== "Low", print(new.data$sex_TraitID), NA )

a$High[!(a$High %in% a$Low)] #6 (previously 4), 208, 2
a$Low[!(a$Low %in% a$High)] #618 (previously 630)


#filter out incomplete entries, and entries with only two data points (will make R2 1.0)
#NOTE: TraitID = 4,630 sex_TraitID = 682,6,208,2,618 (aug 30 two typos need to be fixed, one waiting for email)
new.data<-new.data  %>% 
  filter(!is.na(Treatment)) %>% 
  filter(!is.na(Mean)) %>% 
  filter(PopulationType=="Single") %>% 
  filter(!(sex_TraitID %in% c(682,6,208,2,618))) #these are entries that need to be excluded (only one pred level) 


new.data$Treatment<-as.factor(new.data$Treatment)

#how many traits/studies?
length(unique(new.data$TraitID)) #782 as of Aug 26th (781 without 630)
length(unique(new.data$Study.ID)) #33

#NOTES for sex specific trait IDs ->2 or less levels: 682 (emailed), 208 (typo trait labelled male only sex labelled F),
# 6 (previously trait 4 only one entry), 2 (typo?)

#set up a new function to loop through all data and run anova on Mean trait values for each trait
#note this assumes we have just one covariate of interest: treatment (basically habitat, e.g. high/low predation)
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
    
    aov.mod   <- aov(Mean~Treatment,data=dat[dat$TraitID==levels(trait)[i],])
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


#write results to a csv file for further analyses
#write.table(output.all, file = "TraitR2.csv",row.names=FALSE,col.names=TRUE, sep=",")


