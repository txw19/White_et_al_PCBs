#Analysis - Fish Health- Cortisol
library(ggplot2)
library(jagsUI)
rm(list=ls())

#Read in data
df<-read.table('Fish Health t1_3 Sampling Datasheet 12MAY11_for R.txt',header=T,na.strings='NA')

# make sure all data is imported
head(df)
names(df)

#look at number of rows and column in matrix
dim(df)

#run summary stats
summary(df)

# Make new variable for Sampling Events 
df$time<-ifelse(df$Month==6 & df$Day==2 | df$Day==3,'T1',
ifelse(df$Month==6 & df$Day==29 | df$Day==30,'T2','T3'))
df$time

#Determine if the variable 'Cortisol' is numeric
is.numeric(df$Cortisol)

# Remove NA's from Cortisol column

I<-!is.na(df$Cortisol)
df<-df[I,]
summary(df)
dim(df)
length(unique(df$Cortisol))



#Determine if the variables  are being viewed as factors
is.factor(df$Treatment)
is.factor(df$time)
is.numeric(df$Tank)


#Convert variable 'time' to factors
df$time<-as.factor(df$time)



summary(df)

# Select out columns of interest
df1<-subset(df,select=c('Treatment','time','Cortisol',
'Sex','Tank'))
head(df1)

write.csv(df, 'cort_data_for_matchin.csv')

##############################################
#####Select response variable of interest#####
##############################################

## Obtain summary statistics 
#Activate downloaded library
library(doBy)

sumfun <- function(x, ...){
  c(m=mean(x, ...), sd=sd(x, ...), n=length(x))
}
sum1<-summaryBy(Cortisol~time+Treatment,data=df1,FUN=sumfun,na.rm=T) 
sum1


summary(df1)

#########################
### CREATE BOX PLOT #####
#########################
a<-boxplot(Cortisol~time+Treatment,data=df1)

hist(df1$Cortisol,col='magenta',breaks=30)

#######################
##Log Transform Y Data#
#######################
hist(log(df1$Cortisol),col='magenta',breaks=30)

#Create Object with log transformed ppm values
df1$lcort<-(log(df1$Cortisol))


df1$Tank<-as.factor(df1$Tank)




# Load required package
library(arm)


# Define the model in the BUGS language and write a text file
sink("twowayanovacort1.bug")
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
   cort[i] ~ dnorm(mean[i], tau[trt[i]])    
   mean[i] <- group.mean[trt[i],Time[i]]+ u[tank[i]] #adding tank accounts for tank variability
   } 

#Derived Quantities

tau.tank <- pow(sigma.tank,-2)

} # end model
",fill = TRUE)
sink()


# Bundle data
bugs.data <- list(cort=as.numeric(df1$lcort), Time=as.numeric(df1$time), trt=as.numeric(df1$Treatment),tank=as.numeric(df1$Tank),
n=length(df1$lcort), n.Time=length(unique(df1$time)),n.trt=length(unique(df1$Treatment)),n.tank=length(unique(df1$Tank)))


# Initial values

inits <- function(){list(sigma.tank = rlnorm(1), sigma = rlnorm(4) )}

# , group.mean=rnorm(12)

# Parameters to estimate/keep track of
parameters <- c("group.mean", "sigma", "grp.mean.bt", "sigma.tank", "u")

# MCMC settings
niter <- 10000
nthin <- 1
nburn <- 5000
nchains <- 3



# Do the MCMC stuff calling WinBUGS from R
out2 <- jags(data = bugs.data, inits = inits, parameters.to.save = parameters, 
model.file = "twowayanovacort1.bug", n.chains = nchains, n.thin = nthin, n.iter = niter, 
n.burnin = nburn,parallel = T)


# Summarize the result
print(out2, digits = 3)



