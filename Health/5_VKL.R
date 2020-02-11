#Analysis - Fish Health- Cortisol
library(ggplot2)
library(jagsUI)
library(data.table)
rm(list=ls())

#Read in data
df<-fread('Immune_data_final2.csv')

# make sure all data is imported
head(df)
tail(df)
dim(df)
names(df)

df$TRT <- as.factor(df$TRT)
df$Fish_ID <- as.factor(df$Fish_ID)
df$Time <- as.factor(df$Time)
df$Tank <- as.factor(df$Tank)
#run summary stats
summary(df)


hist(df$VKL)


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

df$y <- log(df$VKL)

# Bundle data
bugs.data <- list(cort=df$y, Time=as.numeric(df$Time), trt=as.numeric(df$TRT),tank=as.numeric(as.factor(as.numeric(df$Tank))),
n=length(df$y), n.Time=length(unique(df$Time)),n.trt=length(unique(df$TRT)),n.tank=length(unique(df$Tank)))


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



