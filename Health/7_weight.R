rm(list=ls())
library(ggplot2)
library(doBy)
library(jagsUI)
library(arm)
library(jagsUI)
# read in data
df<-read.table('Final Gross Examination Data for R.txt',header=T,na.strings='NA')
summary(df)

# Make new variable
df$sample<-ifelse(df$Month==4,'T0',ifelse(df$Month==6 & df$Day==2 | df$Day==3,'T1',
ifelse(df$Month==6 & df$Day==29 | df$Day==30,'T2','T3')))

df$sample <- as.factor(df$sample)

head(df)
summary(df)

sumfun <- function(x, ...){
  c(m=mean(x, ...), sd=sd(x, ...), n=length(x))
}

#Run summary of weight
wgt.sum<-summaryBy(Weight~sample+Treatment,data=df,FUN=sumfun,na.rm=T) 
wgt.sum





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




library(plyr)
df$Treatment <- revalue(df$Treatment, c('Control.Fed'='CF', 'Control.Not.Fed'='CNF','PCB.Fed'= 'PF','PCB.Not.Fed'= 'PNF') )

# Bundle data
bugs.data <- list(ppm=log(df$Weight), Time=as.numeric(df$sample), trt=as.numeric(df$Treatment),tank=as.numeric(as.factor(as.numeric(df$Tank))),
                  n=length(df$Weight), n.Time=length(unique(df$sample)),n.trt=length(unique(df$Treatment)),n.tank=length(unique(df$Tank)))

# Initial values

inits <- function(){list(sigma.tank = runif(1), sigma = runif(4) )}

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


