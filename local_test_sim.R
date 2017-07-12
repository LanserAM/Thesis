#### Survival Simulation Data w/ Low Event Rate and Use of Shrinkage Methods ####
#### by Austin Lanser##############
#### Modified from code by Yuwei Zhu with input by Cindy Chen
#### Last modified: 11/11/2014 ################

rm(list=ls(all=TRUE)) #clear the global environment

library(Hmisc)
library(rms)
library(glmnet)
library(grid)

source("ThesisFunctions.R")

##############################################
#read in data
data<- read.table("Thesis_Simulated_Data.txt", header=TRUE, quote="\"", stringsAsFactors=FALSE)

# input <- read.table("thesis_input.csv", skip=1, header=FALSE, nrows=1)
# i<- as.integer(input[1]) #n1
# j<- as.integer(input[2]) #hazard
# k<- as.integer(input[3]) #exposure rate
# seed <- as.integer(input[4]) 

i<-1
j<-1
k<-1
seed<-859426957

#number of bootstraps
nbt<- 50

# NSIM= number of simulations
NSIM<- 100

## i ## N1v is different numbers of patients considered
N1v <- c(80, 120, 160, 200, 240) 
Kn <- length(N1v)

## j ##log hazards of treatment entertained
Log.Haz<- c(log(1/3), log(1/2), log(1), log(2), log(3))
Kh<- length(Log.Haz)

## k ## Expose.Rate is the proportion of patients with treatment (exposure) 
Expose.Rate <- c(.3 ,.5)
Ktrt<- length(Expose.Rate)

## Simulated survival data will have the following parameters
n1 <- N1v[i]
log.haz<- Log.Haz[j]
rate.expose<- Expose.Rate[k]
tau<- 3
re<- .2
#coefficients when exposure rate of .3 is desired
#phi coefficients are labelled starting w/ 2 so they correspond with coefficients in beta
phi.3<- c(phi2= -1.3, phi3=.44, phi4=.44, phi5=-.88, phi6=-.44, phi7=-.88, phi8=-1.3, phi9=-.44, phi10=-.44, phi11=-.88, phi12=-.44, phi13=-.44, phi14=-.44, phi15=-.44, phi16=-.44, phi17=-.44, phi18=-.44, phi19=-.44, phi20=-.44, phi21=-.44, phi22=-.44, phi23=-.88, phi24=-.44, phi25=-.44, phi26=-.88  ) 

#coefficients when exposure rate of .5 is desired
phi.5<- c(phi2= .9, phi3=-.09, phi4=-.09, phi5=.02, phi6=.02, phi7=.06, phi8=.09, phi9=.02, phi10=.02, phi11=.09, phi12=.02, phi13=.02, phi14=.02, phi15=.02, phi16=.02, phi17=.02, phi18=.02, phi19=.02, phi20=.02, phi21=.02, phi22=.06, phi23=.02, phi24=.02, phi25=.02, phi26=.02  )

#coefficients for event generation
beta1vec<- beta<- c(beta1=log.haz, beta2= -.8, beta3= -8, beta4= -.5, beta5= -.8, beta6= -.8, beta7= -.8, beta8= -.8, beta9= -.3, beta10= -.8, beta11= -.8,  beta12= -.8, beta13= -.8, beta14= -.8, beta15= -.6, beta16= -.8, beta17= -.8, beta18= -.8, beta19= -.7, beta20= -.8, beta21= -.8, beta22= -.8, beta23= -.8, beta24= -.6, beta25= -.8, beta26= -.5, beta27= .4, beta28= .1, beta29= -.3, beta30= .2, beta31= .1)

##make a matrix to keep track of errors
#number of different methods tested
nmethod<- 8
error <- as.data.frame(matrix(0, NSIM, 5+nmethod))
names(error) <- c("N1", "expose.rate", "log.haz", "sim", "Prop", "coxph",  "coxph.PS",  "coxph.Ridge", "coxph.Ridge.nopen", "coxph.Lasso", "coxph.Lasso.nopen", "coxph.unadjust", "coxph.PSHet")
error$N1 <- rep(n1, each=NSIM)
error$expose.rate <- rep(rate.expose, each=NSIM)
error$log.haz <- rep(log.haz, each=NSIM)
error$sim <- c(1:NSIM)

result <- as.data.frame(matrix(NA, NSIM, 4+6*nmethod+8))
names(result) <- c("N1", "expose.rate","log.haz", "sim", "M1.est", "M1.bias", "M1.std", "M1.Cov", "M1.Lower", "M1.Upper","M2.est", "M2.bias", "M2.std", "M2.Cov", "M2.Lower", "M2.Upper", "M3.est", "M3.bias", "M3.mean", "M3.median","M3.std", "M3.Cov", "M3.Lower", "M3.Upper", "M4.est", "M4.bias", "M4.std", "M4.Cov", "M4.Lower", "M4.Upper", "M4.mean", "M4.median","M5.est", "M5.bias", "M5.std", "M5.Cov","M5.Lower", "M5.Upper", "M5.mean", "M5.median","M6.est","M6.bias", "M6.std", "M6.Cov", "M6.Lower", "M6.Upper", "M6.mean", "M6.median", "M7.est", "M7.bias", "M7.std", "M7.Cov", "M7.Lower", "M7.Upper","M8.est", "M8.bias", "M8.std", "M8.Cov", "M8.Lower", "M8.Upper")
result$N1 <- rep(N1v[i], each=NSIM)
result$expose.rate <- rep(Expose.Rate[k], each=NSIM)
result$log.haz <- rep(Log.Haz[j], each=NSIM)
result$sim <- c(1:NSIM)

      set.seed(seed)
      #repeat for number of desired simulations
      estall <- est <- rep(NA, NSIM)
      for (sim in 1:NSIM) {
        
        ########simulate survival data #####
        if(k==1) {
          ret <- surv.sim2(n1, tau, beta1vec=beta, phi1vec=phi.3, rate.event=re, data=data)
          d<- ret$data
          dall <- ret$data0
        }
          if(k==2) d<- surv.sim2(n1, tau, beta1vec=beta, phi1vec=phi.5, rate.event=re, data=data)$data
        
        estall[sim] <- coxph(Surv(ObservedT, Event_Indicator) ~ met.sim +AGE +Female + smoke + tumor_other +
                          lung.sim + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + 
                          cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, 
                        data=dall)$coefficients['met.sim']
        
        est[sim] <- coxph(Surv(ObservedT, Event_Indicator) ~ met.sim +AGE +Female + smoke + tumor_other +
                          lung.sim + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + 
                          cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, 
                          data=d)$coefficients['met.sim']
      }

        print(c(i,j,k,n1,log.haz,rate.expose)); print(Sys.time())
        
        #index for results table
        ind <- sim
        
        ####FOR PROPENSITY SCORE####
        #full list of vars: met.sim +AGE +W +Female + smoke +stage1 +stage2or3 + stage4  + tumor_other +lung.sim + breast + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + cont5 + colorectal.sim + marrow.sim + cervical.sim + hepatic.sim + pancreas.sim + eso_gastric.sim + insulin.sim + uterine.sim + thyroid.sim + oral.sim + ovary.sim + brain.sim
        
        try <- tryCatch(glm(met.sim ~ AGE + Female + smoke + tumor_other + lung.sim + prostate.sim + bladder.sim + BMI + cont1 + cont2 + cont3 + cont4 + cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, data = d, family = "binomial"), error=function(e){1}, warning=function(w){1})
        #I removed: thyroid.sim + eso_gastric.sim + cervical.sim + pancreas.sim + ovary.sim +W +stage1 +stage2or3 +stage4 + hepatic.sim + colorectal + breast
        
        if(!is.list(try)){ 
          d$prop.zscore <- NA
          error[ind,]$Prop <- 1
        }else{
          d$prop.zscore <- as.vector(try$linear.predictors)
        }
        
        ########################################################################
        ## run VE method 1 -- unconditional coxph regression without penalization 
        med0VE <- coxph(Surv(ObservedT, Event_Indicator) ~ met.sim +AGE +Female + smoke + tumor_other +
                                   lung.sim + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + 
                                   cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, data=dall)
        
        med1VE <- coxph(Surv(ObservedT, Event_Indicator) ~ met.sim +AGE +Female + smoke + tumor_other +
                          lung.sim + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + 
                          cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, data=d)
        
        
        med1VE <- tryCatch(coxph(Surv(ObservedT, Event_Indicator) ~ met.sim +AGE +Female + smoke + tumor_other +lung.sim + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, x=TRUE, y=TRUE, data= d), error=function(e){1}, warning=function(w){1} )
        if(!is.list(med1VE)){
          error[ind,"coxph"] <- 1 
        }else{  
          tmp <- summary(med1VE)$coefficients['met.sim',]
          result[ind,"M1.est"] <- tmp['coef']
          result[ind,"M1.bias"] <- tmp['coef'] - log.haz
          result[ind,"M1.std"] <- tmp["se(coef)"]
          conf <- summary(med1VE)$conf.int['met.sim',]
          result[ind,"M1.Lower"] <- conf["lower .95"]
          result[ind,"M1.Upper"] <- conf["upper .95"]
          result[ind,"M1.Cov"] <- ifelse(conf["lower .95"] < exp(log.haz) & conf["upper .95"] > exp(log.haz), 1, 0)
        }
        
        print("M1 Complete"); print(Sys.time())
        
        
        ## run VE method 2 -- unconditional coxph regression with propensity score adjustment ########
        med2VE <- tryCatch(coxph(Surv(ObservedT, Event_Indicator) ~ met.sim + rcs(prop.zscore,4), x=TRUE, y=TRUE, data= d), error=function(e){1}, warning=function(w){1} )
        if(!is.list(med2VE)){
          error[ind,"coxph.PS"] <- 1 
        }else{  
          tmp <- summary(med2VE)$coefficients['met.sim',]
          result[ind,"M2.est"] <- tmp['coef']
          result[ind,"M2.bias"] <- tmp['coef'] - log.haz
          result[ind,"M2.std"] <- tmp["se(coef)"]
          conf<- summary(med2VE)$conf.int['met.sim',]
          result[ind,"M2.Lower"] <- conf["lower .95"]
          result[ind,"M2.Upper"] <- conf["upper .95"]
          result[ind,"M2.Cov"] <- ifelse(conf["lower .95"] < exp(log.haz) & conf["upper .95"] > exp(log.haz), 1, 0)
        }
        print("M2 Complete"); print(Sys.time())
        
        
        ## run VE method 3 -- coxph regression with ridge ########
        x <- as.matrix(model.frame(~ met.sim +AGE +Female + smoke + tumor_other +lung.sim +  prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, data=d))
        tmp <- tryCatch(cv.glmnet(x=x, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=0), error=function(e){1}, warning=function(w){1} )
        med3VE <- tryCatch(glmnet(x=x, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=0, lambda=tmp$lambda.min), error=function(e){1}, warning=function(w){1} )
        
        if(!is.list(med3VE)){
          error[ind,"coxph.Pena"] <- 1 
        }else{        
          result[ind,"M3.est"] <- coef(med3VE)[1]
          result[ind,"M3.bias"] <- coef(med3VE)[1] - log.haz
          
          ## bootstrap to calculate SE ##
          tmprst <- rep(NA, nbt)
          tmpind <- 0
          while (tmpind<nbt) try({
            tempind <- TRUE
            count<- nogood <-0
            while (tempind) {
              tmpdata <- d[sample(1:nrow(d), size=nrow(d), replace=TRUE),]
              tempind <- (sum(table(tmpdata$Event_Indicator, tmpdata$met.sim)>0) < 4)
              if(count>10000) {
                nogood<- error[ind,"coxph.Pena"] <- 1
                break}
              count<- count + 1
            }
            
            if(nogood==1) {
              nogood<-0
              break}
            
            tmpx <- as.matrix(model.frame(~ met.sim +AGE +Female + smoke + tumor_other +lung.sim  + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, data=tmpdata))
            tmp2 <- tryCatch(cv.glmnet(x=x, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=0), error=function(e){1}, warning=function(w){1})
            tmp3 <- tryCatch(glmnet(x=tmpx, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=0, lambda=tmp2$lambda.min), error=function(e){1}, warning=function(w){1})
            
            if(is.list(tmp3)){
              tmpind <- tmpind + 1
              tmprst[tmpind] <- coef(tmp3)[1]
            } 
          }, silent=TRUE)
          
          ## end of bootstrap ##
          
          result[ind,"M3.std"] <- sd(tmprst, na.rm=TRUE)
          result[ind,"M3.Lower"] <- low <- as.numeric(quantile(tmprst, probs=0.025, na.rm=TRUE))
          result[ind,"M3.Upper"] <- upp <- as.numeric(quantile(tmprst, probs=0.975, na.rm=TRUE))
          result[ind,"M3.Cov"] <- ifelse(low < exp(log.haz) & upp > exp(log.haz), 1, 0)
          result[ind, 'M3.mean'] <- mean(tmprst, na.rm=TRUE)
          result[ind, 'M3.median'] <- median(tmprst, na.rm=TRUE)
        }
        print("M3 Complete"); print(Sys.time())
        
        
        ### run VE method 4 -- coxph regression with ridge -- no penalty on metformin status ####### 
        x <- as.matrix(model.frame(~ met.sim +AGE +Female + smoke + tumor_other +lung.sim + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, data=d))
        tmp <- tryCatch(cv.glmnet(x=x, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=0, penalty.factor=c(0,rep(1,ncol(x)-1))), error=function(e){1}, warning=function(w){1})
        med4VE <- tryCatch(glmnet(x=x, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=0, lambda=tmp$lambda.min, penalty.factor=c(0,rep(1,ncol(x)-1))), error=function(e){1}, warning=function(w){1})
        
        if(!is.list(med4VE)){
          error[ind,"coxph.Pena.nopen"] <- 1 
        }else{        
          result[ind,"M4.est"] <- coef(med4VE)[1]
          result[ind,"M4.bias"] <- coef(med4VE)[1] - log.haz
          
          ## bootstrap to calculate SE ##
          tmprst <- rep(NA, nbt)
          tmpind <- 0
          while (tmpind<nbt) try({
            
            tempind <- TRUE
            count<-nogood<-0
            while (tempind) {
              tmpdata <- d[sample(1:nrow(d), size=nrow(d), replace=TRUE),]
              tempind <- (sum(table(tmpdata$Event_Indicator, tmpdata$met.sim)>0) < 4)
              if(count>10000) {
                nogood<- error[ind,"coxph.Pena.nopen"] <- 1
                break}
              count<- count + 1
            }
            
            if(nogood==1) {
              nogood<-0
              break}
            
            tmpx <- as.matrix(model.frame(~ met.sim +AGE +Female + smoke + tumor_other +lung.sim + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, data=tmpdata))
            tmp2 <- tryCatch(cv.glmnet(x=tmpx, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=0, penalty.factor=c(0,rep(1,ncol(tmpx)-1))), error=function(e) {1}, warning=function(w){1})
            tmp3 <- tryCatch(glmnet(x=tmpx, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=0, lambda=tmp2$lambda.min, penalty.factor=c(0,rep(1,ncol(tmpx)-1))), error=function(e){1}, warning=function(w){1})
            
            if(is.list(tmp3)){
              tmpind <- tmpind + 1
              tmprst[tmpind] <- coef(tmp3)[1]
            }
          })
          ## end of bootstrap ##
          
          result[ind,"M4.std"] <- sd(tmprst, na.rm=TRUE)
          result[ind,"M4.Lower"] <- low <- as.numeric(quantile(tmprst, probs=0.025, na.rm=TRUE))
          result[ind,"M4.Upper"] <- upp <- as.numeric(quantile(tmprst, probs=0.975, na.rm=TRUE))
          result[ind,"M4.Cov"] <- ifelse(low < exp(log.haz) & upp > exp(log.haz), 1, 0)
          result[ind, 'M4.mean'] <- mean(tmprst, na.rm=TRUE)
          result[ind, 'M4.median'] <- median(tmprst, na.rm=TRUE)
        }
        print("M4 Complete"); print(Sys.time())
        
        
        ##method 5 -- coxph regression with lasso ########
        x=as.matrix(model.frame(~ met.sim +AGE +Female + smoke + tumor_other +lung.sim + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, data=d))
        #cross-validate for lambda.min: the value of lambda that gives the minimum mean cross-validated error
        tmp <- tryCatch(cv.glmnet(x=x, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=1), error=function(e){1}, warning=function(w){1}) 
        med5VE <- tryCatch(glmnet(x=x, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=1, lambda=tmp$lambda.min), error=function(e){1}, warning=function(w){1})
        
        if(!is.list(med5VE)){
          error[ind,"coxph.Lasso"] <- 1 
        }else{  
          result[ind,"M5.est"] <- coef(med5VE)[1]
          result[ind,"M5.bias"] <- coef(med5VE)[1]- log.haz
          
          ## bootstrap to calculate SE ##
          tmprst <- rep(NA, nbt)
          tmpind <- 0
          while (tmpind<nbt) try({
            tempind <- TRUE
            count<-nogood<-0
            while (tempind) {
              tmpdata <- d[sample(1:nrow(d), size=nrow(d), replace=TRUE),]
              tempind <- (sum(table(tmpdata$Event_Indicator, tmpdata$met.sim)>0) < 4)
              if(count>10000) {
                nogood<- error[ind,"coxph.Lasso"] <- 1
                break}
              count<- count + 1
            }
            
            if(nogood==1) {
              nogood<-0
              break}
            
            tmpx <- as.matrix(model.frame(~ met.sim +AGE +Female + smoke + tumor_other +lung.sim + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, data=tmpdata))
            tmp2 <- tryCatch(cv.glmnet(x=tmpx, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=1), error=function(e) {1}, warning=function(w){1})
            tmp3 <- tryCatch(glmnet(x=tmpx, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=1, lambda=tmp2$lambda.min), error=function(e){1}, warning=function(w){1})
            
            if(is.list(tmp3)){
              tmpind <- tmpind + 1
              tmprst[tmpind] <- coef(tmp3)[1]
            }
          })
          ## end of bootstrap ##
          
          result[ind,"M5.std"] <- sd(tmprst, na.rm=TRUE)
          result[ind,"M5.Lower"] <- low <- as.numeric(quantile(tmprst, probs=0.025, na.rm=TRUE))
          result[ind,"M5.Upper"] <- upp <- as.numeric(quantile(tmprst, probs=0.975, na.rm=TRUE))
          result[ind,"M5.Cov"] <- ifelse(low < exp(log.haz) & upp > exp(log.haz), 1, 0)
          result[ind, 'M5.mean'] <- mean(tmprst, na.rm=TRUE)
          result[ind, 'M5.median'] <- median(tmprst, na.rm=TRUE)
        }
        print("M5 Complete"); print(Sys.time())
        
        
        
        ### run VE method 6 -- coxph regression with Lasso, no penalization on treatment 
        x <- as.matrix(model.frame(~ met.sim +AGE +Female + smoke + tumor_other +lung.sim + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, data=d))
        tmp <- tryCatch(cv.glmnet(x=x, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=1, penalty.factor=c(0,rep(1,ncol(x)-1))), error=function(e){1}, warning=function(w){1})
        med6VE <- tryCatch(glmnet(x=x, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=1, lambda=tmp$lambda.min, penalty.factor=c(0,rep(1,ncol(x)-1))), error=function(e){1}, warning=function(w){1})
  
        
        if(!is.list(med6VE)){
          error[ind,"coxph.Lasso.nopen"] <- 1 
        }else{        
          result[ind,"M6.est"] <- coef(med6VE)[1]
          result[ind,"M6.bias"] <- coef(med6VE)[1] - log.haz
          
          ## bootstrap to calculate SE ##
          tmprst <- rep(NA, nbt)
          tmpind <- 0
          while (tmpind<nbt) try({
            
            tempind <- TRUE
            count<-nogood<-0
            while (tempind) {
              tmpdata <- d[sample(1:nrow(d), size=nrow(d), replace=TRUE),]
              tempind <- (sum(table(tmpdata$Event_Indicator, tmpdata$met.sim)>0) < 4)
              if(count>10000) {
                nogood<- error[ind,"coxph.Lasso.nopen"] <- 1
                break}
              count<- count + 1
            }
            
            if(nogood==1) {
              nogood<-0
              break}
            
            tmpx <- as.matrix(model.frame(~ met.sim +AGE +Female + smoke + tumor_other +lung.sim + prostate.sim +bladder.sim +BMI + cont1 + cont2 + cont3 + cont4 + cont5 + marrow.sim + insulin.sim + uterine.sim + oral.sim + brain.sim, data=tmpdata))
            tmp2 <- tryCatch(cv.glmnet(x=tmpx, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=1, penalty.factor=c(0,rep(1,ncol(tmpx)-1))), error=function(e){1}, warning=function(w){1})
            tmp3 <- tryCatch(glmnet(x=tmpx, y=cbind(time=d$ObservedT, status=d$Event_Indicator), family='cox', alpha=1, lambda=tmp2$lambda.min, penalty.factor=c(0,rep(1,ncol(tmpx)-1))), error=function(e){1}, warning=function(w){1})
            
            if(is.list(tmp3)){
              tmpind <- tmpind + 1
              tmprst[tmpind] <- coef(tmp3)[1]
            }
          }, silent=TRUE)
          ## end of bootstrap ##
          
          result[ind,"M6.std"] <- sd(tmprst, na.rm=TRUE)
          result[ind,"M6.Lower"] <- low <- as.numeric(quantile(tmprst, probs=0.025, na.rm=TRUE))
          result[ind,"M6.Upper"] <- upp <- as.numeric(quantile(tmprst, probs=0.975, na.rm=TRUE))
          result[ind,"M6.Cov"] <- ifelse(low < exp(log.haz) & upp > exp(log.haz), 1, 0)
          result[ind, 'M6.mean'] <- mean(tmprst, na.rm=TRUE)
          result[ind, 'M6.median'] <- median(tmprst, na.rm=TRUE)
        }
        print("M6 Complete"); print(Sys.time())
        
        
        ###### run VE method 7a -- Regular CoxPH unadjusted
        med7VE <- tryCatch(coxph(Surv(ObservedT, Event_Indicator) ~ met.sim, x=TRUE, y=TRUE, data= d), error=function(e){1}, warning=function(w){1})
        if(!is.list(med1VE)){
          error[ind,"coxph.unadjust"] <- 1 
        }else{  
          tmp <- summary(med7VE)$coefficients['met.sim',]
          result[ind,"M7.est"] <- tmp['coef']
          result[ind,"M7.bias"] <- tmp['coef'] - log.haz
          result[ind,"M7.std"] <- tmp["se(coef)"]
          conf <- summary(med1VE)$conf.int['met.sim',]
          result[ind,"M7.Lower"] <- conf["lower .95"]
          result[ind,"M7.Upper"] <- conf["upper .95"]
          result[ind,"M7.Cov"] <- ifelse(conf["lower .95"] < exp(log.haz) & conf["upper .95"] > exp(log.haz), 1, 0)
        }
        
        print("M7 Complete"); print(Sys.time())
        
        
        #         ## run VE method 7 -- unconditional coxph regression with propensity score adjustment after trimming ########
        #         # propensity score (exclude non overlap)
        #         temp1 = d[d$Event_Indicator == 1, ]
        #         temp1Rge = range(temp1$prop.zscore) 
        #         temp0 = d[d$Event_Indicator == 0 & d$prop.zscore >= temp1Rge[1] &  d$prop.zscore <= temp1Rge[2], ]
        #         
        #         if(dim(temp0)[1] >= dim(temp1)[1]){
        #           d01 = rbind(temp0, temp1)
        #           med7VE <- tryCatch(coxph(Surv(ObservedT, Event_Indicator) ~ met.sim + rcs(prop.zscore,4), x=TRUE, y=TRUE, data=d01), error=function(e) 1)
        #           if(!is.list(med7VE)){
        #             error[ind,"coxph.PStrim"] <- 1
        #           }else{  
        #             tmp <- summary(med7VE)$coefficients['met.sim',]
        #             result[ind,"M7.est"] <- tmp['coef']
        #             result[ind,"M7.bias"] <- tmp['coef'] - log.haz
        #             result[ind,"M7.std"] <- tmp["se(coef)"]
        #             result[ind,"M7.Lower"] <- tmp["lower .95"]
        #             result[ind,"M7.Upper"] <- tmp["upper .95"]
        #             result[ind,"M7.Cov"] <- ifelse(tmp["lower 0.95"] < exp(log.haz) & tmp["upper .95"] > exp(log.haz), 1, 0)
        #           }
        #         }else{
        #         cat("    no conv.", "\n")
        #         }
    
        
        ## run VE method 8 -- unconditional coxph regression with propensity score and heterogenecity adjustment########
        med8VE <- tryCatch(coxph(Surv(ObservedT, Event_Indicator) ~ met.sim + rcs(prop.zscore,4) +AGE +smoke +BMI, x=TRUE, y=TRUE, data= d), error=function(e){1}, warning= function(w){1})
        if(!is.list(med8VE)){
          error[ind,"coxph.PSHet"] <- 1 
        }else{  
          tmp <- summary(med8VE)$coefficients['met.sim',]
          result[ind,"M8.est"] <- tmp['coef']
          result[ind,"M8.bias"] <- tmp['coef'] - log.haz
          result[ind,"M8.std"] <- tmp["se(coef)"]
          conf<- summary(med8VE)$conf.int['met.sim',]
          result[ind,"M8.Lower"] <- conf["lower .95"]
          result[ind,"M8.Upper"] <- conf["upper .95"]
          result[ind,"M8.Cov"] <- ifelse(conf["lower .95"] < exp(log.haz) & conf["upper .95"] > exp(log.haz), 1, 0)
        }
        ############################################################################################
        #write.table(result, file="test_result.txt", row.names=FALSE, col.names=FALSE)
        #write.table(error, file="test_error.txt", row.names=FALSE, col.names=FALSE)
        cat("k=", k, "i=", i, "j=", j, "sim=", sim,  "  "); print(Sys.time())
      } #end of for/sim 
    

#write.table(result, file="test_result.txt", row.names=FALSE, col.names=FALSE)
#write.table(error, file="test_error.txt", row.names=FALSE, col.names=FALSE)