#Final Project! QME 
library(ggplot2)
rm(list=ls())

#Read in data
df<-read.table('Total PCBs 2FEB11_HW2_29APR11.txt',header=T,na.strings='NA')

# make sure all data is imported
head(df, 20)
names(df)

#look at number of rows and column in matrix
dim(df)

#run summary stats
summary(df)

#Create a new variable for ppm
df$ppm<-df$Result/1000000

#Determine if the variable 'Result' is numeric
is.numeric(df$Result)

#Determine if the variables 'Client.ID', 'Time', 'Units', 'Sample', and 'Tissue' are being viewed as factors
is.factor(df$Client.ID)
is.factor(df$Time)
is.factor(df$Units)
is.factor(df$sample)
is.factor(df$tissue)

#Convert variables 'Sample' and 'Tissue 'to factors
df$sample<-as.factor(df$sample)
df$tissue<-as.factor(df$tissue)


summary(df)

# Select out columns of interest
df2<-subset(df,select=c('Client.ID','Time','Compound',
'Result','Units','ppm'))
head(df2)


# Look at all unique compound names
 unique(df2$Compound)


##############################################
#####Select response variable of interest#####
##############################################

df3<-subset(df2,Compound=='Total.PCBs')
name<-'Total PCB (mg/kg)'
summary(df3)
head(df3)

# Extract the first 3 letters from a character string
#load library
library(stringr)
# Create object that contains what to extract using the "str_extract" function below
regexp <- "([[:alpha:]]{3})"

# Call a new variable 'sample' that contains the 1st 3 characters of Client.ID
df3$sample<-str_extract(df3$Client.ID,regexp)
head(df3)

# Extract the first 2 letters from a character string
library(stringr)
# Create object that contains what to extract using the "str_extract" function below
regexp <- "([[:alpha:]]{2})"

# Call a new variable 'trt' that contains the 1st 2 characters of Sample
df3$trt<-str_extract(df3$sample,regexp)
head(df3)

# Extract tissue type and make new column for trt
regexp2 <- "([[:lower:]])"

# Call a new variable 'tissue' that contains the lower case value located in 'sample' create above
df3$tissue<-str_extract(df3$sample,regexp2)
head(df3)

# Obtain Tank IDs and place into new column
df3$tank<-substr(df3$Client.ID,5,7)
head(df3)
unique(df3$tank)

# Select a tissue type to examine
df4<-subset(df3,tissue=='w')
head(df4)
dim(df4)
unique(df4$tank)

## Obtain summary statistics 
#Activate downloaded library
library(doBy)

sumfun <- function(x, ...){
  c(m=mean(x, ...), sd=sd(x, ...), n=length(x))
}
sum1<-summaryBy(ppm~Time+trt,data=df4,FUN=sumfun,na.rm=T) 
sum1


#############
## Effects ##
#############

#Make sure effects are viewed correctly for analysis


df4$tank<-as.factor(df4$tank)
df4$sample<-as.factor(df4$sample)
df4$trt<-as.factor(df4$trt)

summary(df4)

#########################
### CREATE BOX PLOT #####
#########################
a<-boxplot(ppm~Time+trt,data=df4)


table(df4$sample,df4$tank)
summary(df4)

hist(df4$ppm,col='magenta',breaks=30)

#######################
##Log Transform Y Data#
#######################
hist(log(df4$ppm),col='magenta',breaks=30)

#Create Object with log transformed ppm values
df4$l.ppm<-(log(df4$ppm))


###########################
## Means Parameterization##
##		of 	       ##
##	Frequentist Model  ##
##				 ##
###########################

# Effects parameterization
# lm.ep<-lm(l.ppm~ Time + trt + Time:trt, data=df4)
# lm.ep
# summary(lm.ep)
# # plot(lm.ep)
# 
# #Means parameterization
# lm.mp<-lm(l.ppm~  Time*trt-1-trt-Time, data=df4)
# lm.mp
# summary(lm.mp)
# plot(lm.mp)

##########################################
## 2 Way ANOVA with Back Transformation ##
##			and                   ##
##    Random Tank Effect in Winbugs     ##
##########################################

# Load required package
library(arm)
library(jagsUI)

# Define the model in the BUGS language and write a text file
sink("twowayanovaPCB.bug")
cat("
model {

for(i in 1:n.tank){
	u[i] ~ dnorm(0, tau.tank) #add tank as a random effect
}

# Priors
for(i in 1:n.trt){
	for(j in 1:n.Time){
	group.mean[i,j]~dnorm(0,0.0001) # normal priors for each trt-time mean
	grp.mean.bt[i,j] <- exp(group.mean[i,j])        # back-transfrom from log-scale to natural scale
	}          
	
}

# Priors for random tank effect

sigma.tank ~ dunif(0,100) # SD hyperparameter for random intercepts (prior of variance among intercepts)
 
for(i in 1:n.trt){
	sigma[i] ~ dunif(0, 100)
	tau[i]<- pow(sigma[i],-2)
}

# Likelihood: 
for (i in 1:n){
   ppm[i] ~ dnorm(mean[i], tau[trt[i]])    
   mean[i] <- group.mean[trt[i],Time[i]]+ u[tank[i]] #adding tank accounts for tank variability
   } 

#Derived Quantities

tau.tank <- pow(sigma.tank,-2)

} # end model
",fill = TRUE)
sink()


# Bundle data
bugs.data <- list(ppm=as.numeric(df4$l.ppm), Time=as.numeric(df4$Time), trt=as.numeric(df4$trt),tank=as.numeric(df4$tank),
n=length(df4$l.ppm), n.Time=length(unique(df4$Time)),n.trt=length(unique(df4$trt)),n.tank=length(unique(df4$tank)))

# Initial values

inits <- function(){list(sigma.tank = rlnorm(1), sigma = rlnorm(4) )}

# , group.mean=rnorm(16)

# Parameters to estimate/keep track of
parameters <- c("group.mean", "sigma", "grp.mean.bt", "sigma.tank", "u")

# MCMC settings
niter <- 10000
nthin <- 1
nburn <- 5000
nchains <- 3


# bugs.dir <-'C:\\WinBUGS14'



# Do the MCMC stuff calling WinBUGS from R
out2 <- jags(data = bugs.data, inits = inits, parameters.to.save = parameters, 
model.file = "twowayanovaPCB.bug", n.chains = nchains, n.thin = nthin, n.iter = niter, 
n.burnin = nburn, parallel = T)


# Summarize the result
print(out2, digits = 3)




