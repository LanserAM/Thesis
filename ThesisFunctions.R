#### Setup for Thesis Simulation  ####
#### by Austin Lanser##############
#### Modified from code by Yuwei Zhu with input by Cindy Chen
#### Last modified: 11/11/2014 ################

library(Hmisc)
library(rms)
library(glmnet)
library(grid)

#############################

makeContinuousList = function(x, minVal, maxVal, precision = 3, distribution = rnorm, seed = NULL, ...)
{
  if(!is.null(seed)){
    set.seed(seed)
  }
  dist = distribution(2* x, ...)
  minD = min(dist)
  maxD = max(dist)
  vSpan = maxVal-minVal
  val = round(vSpan/(maxD-minD)*(dist-minD), digits = precision)
  val = minVal+val
  sample(val,x, replace = TRUE)
}

##############

##sample from the metformin data to get covariates, without replacement to maintain correlation structure
surv.sim2<- function(n, tau, beta1vec, phi1vec, rate.event, data){
  # This simulates survival data
  # n: number of subjects
  # tau: the length of study/cutoff for follow-up 
  # beta1vec: the vector of regression coefficients
  # phi1vec: the vector of regression coefficients used to generate exposure
  # rate.event: the ratio of events to non-events in the final sample
  
  #example parameters
  #tau <- 3
  #beta1vec <- c(-0.5, 1, 0)
  #phi1vec <- c(0, -0.5, 1)
  #rate.event<- .05
  
  #only use one time point per subject (not longitudinal)
  data<- data[!duplicated(data$RUID), ]   
  
  #create replacement binary vars for covariates with under 9% to avoid singularity issues later
  phi<- .25* c(phi2= -1.3, phi3=.44, phi4=.44, phi5=-.88, phi6=-.44, phi7=-.88, phi8=-1.3, phi9=-.44, phi10=-.44, phi11=-.88, phi12=-.44, phi13=-.44, phi14=-.44, phi15=-.44, phi16=-.44, phi17=-.44, phi18=-.44, phi19=-.44, phi20=-.44, phi21=-.44, phi22=-.44, phi23=-.88, phi24=-.44, phi25=-.44, phi26=-.88  )
  
  colorectal.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:7],10,phi[9:25]))))
  lung.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:9],10,phi[11:25]))))
  breast.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:10],10,phi[12:25]))))
  prostate.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:11],10,phi[13:25]))))
  marrow.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:12],10,phi[14:25]))))
  uterine.sim <- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:13],10,phi[15:25]))))
  cervical.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:14],10,phi[16:25]))))
  bladder.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:15],10,phi[17:25]))))
  brain.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:16],10,phi[18:25]))))
  thyroid.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:17],10,phi[19:25]))))
  pancreas.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:18],10,phi[20:25]))))
  ovary.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:19],10,phi[21:25]))))
  hepatic.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:20],10,phi[22:25]))))
  eso_gastric.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:21],10,phi[23:25]))))
  oral.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:22],10,phi[24:25]))))
  insulin.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -c(phi[1:23],10,phi[25]))))
  
  
  #create additional continuous random variables
  cont1 <- round(makeContinuousList(nrow(data), minVal=0, maxVal=100, precision = 3, distribution = rnorm, sd= 5, mean =35),0)
  cont1 <- (cont1-mean(cont1))/sd(cont1)
  
  cont2 <- round(makeContinuousList(nrow(data), minVal=0, maxVal=365, precision = 3, distribution = rnorm, sd= 10, mean =100),0)
  cont2 <- (cont2-mean(cont2))/sd(cont2)
  
  cont3 <- makeContinuousList(nrow(data), minVal=1.5, maxVal=3, precision = 3, distribution = rnorm, sd= .5, mean =1.9)
  cont3 <- (cont3-mean(cont3))/sd(cont3)
  
  cont4 <- round(makeContinuousList(nrow(data), minVal=70, maxVal=500, precision = 3, distribution = rnorm, sd= 30, mean =140),0)
  cont4 <- (cont4-mean(cont4))/sd(cont4)
  
  cont5 <- round(makeContinuousList(nrow(data), minVal=-10, maxVal=10, precision = 3, distribution = rnorm, sd= 1, mean= 0),0)
  cont5 <- (cont5-mean(cont5))/sd(cont5)
  
  #simulate exposure (metformin use)
  met.sim<- rbinom(nrow(data), size=1, prob= 1/(1+exp(as.matrix(data[,2:26]) %*% -phi1vec)))
  mean(met.sim)
  
  data2 <- cbind(met.sim, data[,2:8], colorectal.sim, tumor_other=data[,10], lung.sim, breast=data[,12], prostate.sim, marrow.sim, 
                 uterine.sim, cervical.sim, bladder.sim, brain.sim, thyroid.sim, pancreas.sim, ovary.sim, hepatic.sim, eso_gastric.sim, 
                 oral.sim, insulin.sim, 'BMI'=data[,26], cont1, cont2, cont3, cont4, cont5)
  
  # Te is the true underlying event time
  # baseline hazard is exp(t), cumulative baseline hazard e^(t)-1;
  #Te <- log(-log(runif(nrow(data2),0,1))*exp(-(as.matrix(data2))%*%beta1vec)+1)
  # Another example: baseline hazard 2t, cumulative baseline hazard (t^2)
  #Te <- sqrt(-log(runif(nrow(data2),0,1))*exp(-(as.matrix(data2))%*%beta1vec))
  
  #parametric example using weibull distribution
  Te <- rweibull(nrow(data2), shape=10, scale=1.24*exp(-(as.matrix(data2))%*%(beta1vec/10)))
  summary(Te)
  #min(Te)
  #hist(Te, breaks=1000)
  # generate censoring time, unif(0,6) and then truncate by tau
  Cens <- 7*runif(nrow(data2))
  Cens[Cens>tau] <- tau
  
  # Observed data: Delta is the event indicator and W is observed time
  DeltaE <- ifelse(Te <= Cens, 1, 0)
  mean(DeltaE)
  W <- ifelse(DeltaE==1, Te, Cens)
  
  #make a data frame to work with
  Xdat <- as.data.frame(cbind(data$RUID , round(W, digits=4), DeltaE, data2))
  names(Xdat)[1]<- "RUID"
  names(Xdat)[2]<- "ObservedT"
  names(Xdat)[3]<- "Event_Indicator"
  
  #draw sample id for events
  sampid1<- sample(Xdat$RUID[Xdat$Event_Indicator==1], size=n*rate.event, replace=FALSE)
  
  #draw sample data from using the sampled ids
  sampledata1<- subset(Xdat, Xdat$RUID %in% sampid1)
  
  #draw sample id for non-events
  sampid2<- sample(Xdat$RUID[Xdat$Event_Indicator==0], size=n*(1-rate.event), replace=FALSE)
  
  #draw sample data using sampled ids
  sampledata2<- subset(Xdat, Xdat$RUID %in% sampid2)
  
  Xdat2<- rbind(sampledata1, sampledata2)
  
  #percent of events
  p<- mean(DeltaE)
  a<- mean(met.sim)
  a2<- mean(Xdat2$met.sim)
  return(list(data0=Xdat, data=Xdat2, prop.events=p, true.events=Te, prop.exposure=a, sample.exposure=a2))
}