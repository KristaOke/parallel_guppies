library(stats)

#setwd("<YOURPATHWAY>")
wd <- getwd()
new.data<- read.csv(paste(wd,'/data/YOUR_DATA.csv',sep=""), header=TRUE, sep=",") %>% #let's update with the real file name
  rename( Mean = MeanValue, Treatment = Predation)

#my data was read in as a csv file with columns: 

# "Study.ID"             "Species"              "Contrast"             "Trait.ID"             "Trait"                "Trait.type"           "Standardized."       
# "Standardization.type" "Paired."              "Study.type"           "Data.point"           "Treatment"            "Population"           "Mean"                
# "N" 



#filter out incomplete entries, and entries with only two data points (will make R2 1.0)
new.data<-new.data  %>% 
  filter(!is.na(Treatment)) %>% 
  filter(!is.na(Mean)) %>% 
  filter(PopulationType=="Single") %>% 
  filter(!(TraitID %in% c(277, 282, 283, 284, 285, 326))) #277 only has low, 282-285 are mixed but labelled single, 326 only has one for each but is single?

new.data$Treatment<-as.factor(new.data$Treatment)

#how many traits/studies?
length(unique(new.data$TraitID))
length(unique(new.data$Study.ID))


#set up a new function to loop through all data and run anova on Mean trait values for each trait
#note this assumes we have just one covariate of interest: treatment (basically habitat, e.g. high/low predation)
#if we want more covariates we'll need to rethink this a bit

#START FUNC

dat<-FALSE 

i <- 1

anova_loop <- function(dat=dat){
  
  trait <- as.factor(dat$TraitID)
  
  output <- matrix(ncol=7,nrow=length(levels(trait)))
  colnames(output) <- c("R^2","adj.R^2","partial.Eta","F","p","traitID", "justr")
  
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
    traitID <- levels(trait)[i]
    justr <- sqrt(anova.mod[1,4]/(anova.mod[1,4]+anova.mod[2,1]))
    
    output[i,] <- c(R,R_adj,partial_eta,F_value,p_value, traitID, justr) 
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



