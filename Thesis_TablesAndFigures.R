library("lattice")

# error1 <- read.table("~/Thesis/error1.txt", quote="\"", stringsAsFactors=FALSE)
# names(error) <- c("N1", "expose.rate", "log.haz", "sim", "Prop", "coxph",  "coxph.PS",  "coxph.Ridge", "coxph.Ridge.nopen", "coxph.Lasso", "coxph.Lasso.nopen", "coxph.PStrim", "coxph.PSHet")

###TEST###
total <- read.table("~/Thesis/testresult.txt", quote="\"", stringsAsFactors=FALSE)
names(total) <- c("N1", "expose.rate","log.haz", "sim", "M1.est", "M1.bias", "M1.std", "M1.Cov", "M1.Lower", "M1.Upper","M2.est", "M2.bias", "M2.std", "M2.Cov", "M2.Lower", "M2.Upper", "M3.est", "M3.bias", "M3.mean", "M3.median","M3.std", "M3.Cov", "M3.Lower", "M3.Upper", "M4.est", "M4.bias", "M4.std", "M4.Cov", "M4.Lower", "M4.Upper", "M4.mean", "M4.median","M5.est", "M5.bias", "M5.std", "M5.Cov","M5.Lower", "M5.Upper", "M5.mean", "M5.median","M6.est","M6.bias", "M6.std", "M6.Cov", "M6.Lower", "M6.Upper", "M6.mean", "M6.median", "M7.est", "M7.bias", "M7.std", "M7.Cov", "M7.Lower", "M7.Upper","M8.est", "M8.bias", "M8.std", "M8.Cov", "M8.Lower", "M8.Upper")
############################################################################

result1 <- read.table("~/Thesis/result1.txt", quote="\"", stringsAsFactors=FALSE)

result2 <- read.table("~/Thesis/result2.txt", quote="\"", stringsAsFactors=FALSE)

result3 <- read.table("~/Thesis/result3.txt", quote="\"", stringsAsFactors=FALSE)

result4 <- read.table("~/Thesis/result4.txt", quote="\"", stringsAsFactors=FALSE)

result5 <- read.table("~/Thesis/result5.txt", quote="\"", stringsAsFactors=FALSE)

result6 <- read.table("~/Thesis/result6.txt", quote="\"", stringsAsFactors=FALSE)

result_extra1 <-  read.table("~/Thesis/result_extra1.txt", quote="\"", stringsAsFactors=FALSE)

result_extra2 <- read.table("~/Thesis/result_extra2.txt", quote="\"", stringsAsFactors=FALSE)

total<-rbind(result1, result2, result3, result4, result5, result6)
names(total) <- c("N1", "expose.rate","log.haz", "sim", "M1.est", "M1.bias", "M1.std", "M1.Cov", "M1.Lower", "M1.Upper","M2.est", "M2.bias", "M2.std", "M2.Cov", "M2.Lower", "M2.Upper", "M3.est", "M3.bias", "M3.mean", "M3.median","M3.std", "M3.Cov", "M3.Lower", "M3.Upper", "M4.est", "M4.bias", "M4.std", "M4.Cov", "M4.Lower", "M4.Upper", "M4.mean", "M4.median","M5.est", "M5.bias", "M5.std", "M5.Cov","M5.Lower", "M5.Upper", "M5.mean", "M5.median","M6.est","M6.bias", "M6.std", "M6.Cov", "M6.Lower", "M6.Upper", "M6.mean", "M6.median", "M7.est", "M7.bias", "M7.std", "M7.Cov", "M7.Lower", "M7.Upper","M8.est", "M8.bias", "M8.std", "M8.Cov", "M8.Lower", "M8.Upper")



################## Missing ###############
#assign each row 1 if there are any na's or 0 if complete (not counting method 7)
total$AnyNA<- ifelse(!complete.cases(total),1,0)
nasubset<- total[,c("N1", "expose.rate", "log.haz", "AnyNA")]
nacount<- aggregate(nasubset$AnyNA, list(nasubset$N1, nasubset$expose.rate, nasubset$log.haz), sum)
names(nacount)<- c("N1", "expose.rate", "log.haz", "SumNA")
nacount$log.haz<- factor(nacount$log.haz, levels=c(log(1/3),log(1/2),log(1),log(2),log(3)), labels=c("Log Hazard=log(1/3)", "Log Hazard=log(1/2)", "Log Hazard=log(1)", "Log Hazard= log(2)", "Log Hazard= log(3)"))
nacount$N1<- factor(nacount$N1, levels=c(80,120,160,200,240))

fig.na.3<- barchart(SumNA~N1|log.haz,
                    horizontal=FALSE,
                    data=subset(nacount, nacount$expose.rate==.3),
                    main="Figure: Count NAs (Exposure=.3)", 
                    ylab="Number of NAs", xlab="Sample Size")
fig.na.3

fig.na.5<- barchart(SumNA/max(total$sim)~N1|log.haz,
                    horizontal=FALSE,
                    data=subset(nacount, nacount$expose.rate==.5),
                    main="Figure: Ratio NAs (Exposure=.5)",
                    ylab="Ratio NAs to Total", xlab="Sample Size")
fig.na.5

####### Data setup to do NA figures separated by method
total$case.per.par<- .2*total$N1/18
total$haz<- exp(total$log.haz)
#name methods order ( 1, 7, 2, 8, 3, 4, 5, 6)
mtds<- c("Regular CoxPH", "Regular Unadjusted CoxPH", "PS", "PS Hetero", "Ridge", "Ridge(No Penalty on Exposure)", "LASSO", "LASSO(No Penalty on Exposure)")

total2 <- data.frame(N=rep(total$N1,8), haz=rep(total$haz,8), log.haz=rep(total$log.haz,8), expose.rate=rep(total$expose.rate,8), bias=as.numeric(unlist(total[c('M1.bias','M7.bias','M2.bias','M8.bias', 'M3.bias','M4.bias','M5.bias', 'M6.bias')])), Lower=as.numeric(unlist(total[c('M1.Lower','M7.Lower','M2.Lower', 'M8.Lower', 'M3.Lower','M4.Lower','M5.Lower','M6.Lower')])), Upper=as.numeric(unlist(total[c('M1.Upper','M7.Upper','M2.Upper', 'M8.Upper', 'M3.Upper','M4.Upper','M5.Upper','M6.Upper')])), Coverage=as.numeric(unlist(total[c('M1.Cov','M7.Cov','M2.Cov','M8.Cov','M3.Cov','M4.Cov','M5.Cov','M6.Cov')])), std=as.numeric(unlist(total[c('M1.std','M7.std', 'M2.std', 'M8.std', 'M3.std','M4.std','M5.std','M6.std')])), est=as.numeric(unlist(total[c('M1.est','M7.est' ,'M2.est', 'M8.est', 'M3.est','M4.est','M5.est','M6.est')])), method=rep(mtds, each=nrow(total)))

total2$AnyNA<- ifelse(!complete.cases(total2),1,0)
total3 <- with(total2, aggregate(AnyNA, by=list(N, expose.rate, log.haz, method), FUN=function(x) {
  sum(x, na.rm=TRUE)
}))
colnames(total3) <- c('N', 'expose.rate', 'log.haz', 'method', 'SumNA')

total3$log.haz<- factor(total3$log.haz, levels=c(log(1/3),log(1/2),log(1),log(2),log(3)), labels=c("Log Hazard=log(1/3)", "Log Hazard=log(1/2)", "Log Hazard=log(1)", "Log Hazard= log(2)", "Log Hazard= log(3)"))
total3$N<- factor(total3$N, levels=c(80,120,160,200,240))

fig.na.3.Regular<- barchart(SumNA~N|log.haz,
                    horizontal=FALSE,
                    data=subset(total3, total3$expose.rate==.3 & total3$method=="Regular CoxPH"),
                    index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                    main="Figure: Count NAs for Regular CoxPH (Exposure=.3)", 
                    ylab="Number of NAs", xlab="Sample Size")
fig.na.3.Regular

fig.na.5.Regular<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.5 & total3$method=="Regular CoxPH"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for Regular CoxPH (Exposure=.5)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.5.Regular


#### PS
fig.na.3.PS<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.3 & total3$method=="PS"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for PS (Exposure=.3)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.3.PS

fig.na.5.PS<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.5 & total3$method=="PS"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for PS (Exposure=.5)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.5.PS

#### PS Hetero
fig.na.3.PSH<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.3 & total3$method=="PS Hetero"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for PS Hetero (Exposure=.3)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.3.PSH

fig.na.5.PSH<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.5 & total3$method=="PS Hetero"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for PS Hetero (Exposure=.5)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.5.PSH

#### Ridge
fig.na.3.Ridge<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.3 & total3$method=="Ridge"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for Ridge (Exposure=.3)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.3.Ridge

fig.na.5.Ridge<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.5 & total3$method=="Ridge"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for Ridge (Exposure=.5)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.5.Ridge

#### Ridge No Pen
fig.na.3.RidgeNP<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.3 & total3$method=="Ridge(No Penalty on Exposure)"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for Ridge No Pen (Exposure=.3)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.3.RidgeNP

fig.na.5.RidgeNP<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.5 & total3$method=="Ridge(No Penalty on Exposure)"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for Ridge No Pen (Exposure=.5)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.5.RidgeNP

#### LASSO
fig.na.3.LASSO<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.3 & total3$method=="LASSO"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for LASSO (Exposure=.3)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.3.LASSO

fig.na.5.LASSO<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.5 & total3$method=="LASSO"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for LASSO (Exposure=.5)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.5.LASSO

#### LASSONP
fig.na.3.LASSONP<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.3 & total3$method=="LASSO(No Penalty on Exposure"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for LASSO No Pen (Exposure=.3)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.3.LASSONP

fig.na.5.LASSONP<- barchart(SumNA~N|log.haz,
                            horizontal=FALSE,
                            data=subset(total3, total3$expose.rate==.5 & total3$method=="LASSO(No Penalty on Exposure"),
                            index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                            main="Figure: Count NAs for LASSO No Pen (Exposure=.5)", 
                            ylab="Number of NAs", xlab="Sample Size")
fig.na.5.LASSONP


################# Subset of complete simulations (no NA) ################
total<- total[complete.cases(total)]
################# Bias of Estimates ######
bias<- total[,c("N1", "expose.rate", "log.haz", "sim", "M1.bias", "M2.bias", "M3.bias", "M4.bias", "M5.bias", "M6.bias", "M7.bias", "M8.bias")]
bias.mean<- aggregate(bias[,-c(1,2,3,4)], list(bias$N1, bias$expose.rate, bias$log.haz), mean)
names(bias.mean)[c(1,2,3)]<- c("N1","expose.rate","log.haz")
case.per.par<- .2*bias.mean$N1/18
bias.mean<- cbind(case.per.par, bias.mean)
bias.mean$log.haz<- factor(bias.mean$log.haz, levels=c(log(1/3),log(1/2),log(1),log(2),log(3)), labels=c("Log Hazard=log(1/3)", "Log Hazard=log(1/2)", "Log Hazard=log(1)", "Log Hazard= log(2)", "Log Hazard= log(3)"))
bias.mean<-bias.mean[with(bias.mean, order(case.per.par)), ]

###key###
mykey <- list(space="top", text=list(c("Regular CoxPH","Regular Unadjusted CoxPH", "PS", "PS Hetero", "Ridge", "Ridge(No Penalty on Exposure)", "LASSO", "LASSO(No Penalty on Exposure)")), lines=list(lty=c(1,2,1,2,1,2,1,2), col=c('black', 'black', 'orange', 'orange', 'blue', 'blue', 'green', 'green')), columns=4)

#line graph for combination of log hazard and exposure rate 
pdf("Bias.pdf", height=9, width=9)
fig.bias<- xyplot(M1.bias+M7.bias+M2.bias+M8.bias+M3.bias+M4.bias+M5.bias+M6.bias~case.per.par|log.haz,
       data=bias.mean,
       type= c("l"),
       index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
       par.settings=list(superpose.line=list(col=c('black', 'black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), lty=c(1,2,1,2,1,2,1,2))),
                         #superpose.symbol=list(col=c('black', 'black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'))),
       abline=list(h=0, col='grey'),
       key=mykey,
       main="Figure1: Bias of Estimates", 
       ylab="Bias", xlab="# Cases per Parameter")
fig.bias
dev.off()
################# bias graphs separated by exposure rate
pdf("Bias_3.pdf", height=9, width=9)
fig.bias.3<- xyplot(M1.bias+M7.bias+M2.bias+M8.bias+M3.bias+M4.bias+M5.bias+M6.bias~case.per.par|log.haz,
                  data=subset(bias.mean, bias.mean$expose.rate==.3),
                  index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                  type= c("l"),
                  par.settings=list(superpose.line=list(col=c('black', 'black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), lty=c(1,2,1,2,1,2,1,2))),
                  #superpose.symbol=list(col=c('black', 'black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'))),
                  abline=list(h=0, col='grey'),
                  key=mykey,
                  main="Figure: Bias of Estimates (Exposure=.3)", 
                  ylab="Bias", xlab="# Cases per Parameter")
fig.bias.3
dev.off()

pdf("Bias_5.pdf", height=9, width=9)
fig.bias.5<- xyplot(M1.bias+M7.bias+M2.bias+M8.bias+M3.bias+M4.bias+M5.bias+M6.bias~case.per.par|log.haz,
                    data=subset(bias.mean, bias.mean$expose.rate==.5),
                    index.cond=list(c(1,2,3,5,4)),#this provides the order of the panels
                    type= c("l"),
                    par.settings=list(superpose.line=list(col=c('black', 'black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), lty=c(1,2,1,2,1,2,1,2))),
                    #superpose.symbol=list(col=c('black', 'black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'))),
                    abline=list(h=0, col='grey'),
                    key=mykey,
                    main="Figure: Bias of Estimates (Exposure=.5)", 
                    ylab="Bias", xlab="# Cases per Parameter")
fig.bias.5
dev.off()

#########################Data set up################################################
total$case.per.par<- .2*total$N1/18
total$haz<- exp(total$log.haz)
#name methods order ( 1, 7, 2, 8, 3, 4, 5, 6)
mtds<- c("Regular CoxPH", "Regular Unadjusted Coxph", "PS", "PS Hetero", "Ridge", "Ridge(No Penalty on Exposure)", "LASSO", "LASSO(No Penalty on Exposure)")

total2 <- data.frame(N=rep(total$N1,8), haz=rep(total$haz,8), log.haz=rep(total$log.haz,8), expose.rate=rep(total$expose.rate,8), bias=as.numeric(unlist(total[c('M1.bias','M7.bias','M2.bias','M8.bias', 'M3.bias','M4.bias','M5.bias', 'M6.bias')])), Lower=as.numeric(unlist(total[c('M1.Lower','M7.Lower','M2.Lower', 'M8.Lower', 'M3.Lower','M4.Lower','M5.Lower','M6.Lower')])), Upper=as.numeric(unlist(total[c('M1.Upper','M7.Upper','M2.Upper', 'M8.Upper', 'M3.Upper','M4.Upper','M5.Upper','M6.Upper')])), Coverage=as.numeric(unlist(total[c('M1.Cov','M7.Cov','M2.Cov','M8.Cov','M3.Cov','M4.Cov','M5.Cov','M6.Cov')])), std=as.numeric(unlist(total[c('M1.std','M7.std', 'M2.std', 'M8.std', 'M3.std','M4.std','M5.std','M6.std')])), est=as.numeric(unlist(total[c('M1.est','M7.est', 'M2.est', 'M8.est', 'M3.est','M4.est','M5.est','M6.est')])), method=rep(mtds, each=nrow(total)))

#Lower2 and Upper2 are 1 if they are on the correct side of the hazard
total2$Lower2 <- ifelse(total2$Lower <= total2$haz, 1, 0)
total2$Upper2 <- ifelse(total2$Upper >= total2$haz, 1, 0)
total2$case.per.par<- .2*total$N/18

total3 <- with(total2, aggregate(bias, by=list(N, expose.rate, log.haz, method), FUN=mean, na.rm=TRUE))
colnames(total3) <- c('N', 'expose.rate', 'log.haz', 'method', 'bias')

total3$mse <- with(total2, aggregate(bias, by=list(N, expose.rate, log.haz, method), FUN=function(x) {
  mean(x^2, na.rm=TRUE)
}))$x

total3$mae <- with(total2, aggregate(bias, by=list(N, expose.rate, log.haz, method), FUN=function(x) {
  median(abs(x), na.rm=TRUE)
}))$x

total3$Coverage <- with(total2, aggregate(Coverage, by=list(N, expose.rate, log.haz, method), FUN=mean, na.rm=TRUE))$x

total3$meanae <- with(total2, aggregate(bias, by=list(N, expose.rate, log.haz, method), FUN=function(x) {
  mean(abs(x), na.rm=TRUE)
}))$x

total3$lower.cov <- with(total2, aggregate(Lower2, by=list(N, expose.rate, log.haz, method), FUN=function(x) {
  mean(x, na.rm=TRUE)
}))$x

total3$upper.cov <- with(total2, aggregate(Upper2, by=list(N, expose.rate, log.haz, method), FUN=function(x) {
  mean(x, na.rm=TRUE)
}))$x

# empirical se #
total3$ese <- with(total2, aggregate(est, by=list(N, expose.rate, log.haz, method), FUN=sd, na.rm=TRUE))$x

total3$method <- factor(total3$method, levels=c("Regular CoxPH", "Regular Unadjusted CoxPH", "PS", "PS Hetero", "Ridge", "Ridge(No Penalty on Exposure)", "LASSO", "LASSO(No Penalty on Exposure)"))
total3$log.haz.lab <- paste('Log Hazard =', format(round(total3$log.haz,2), nsmall=2))

total3$case.per.par<- .2*total3$N/18

total3 <- total3[order(total3$case.per.par),]


#########################Empirical Standard Error of Estimates
pdf("ESE.pdf", height=9, width=9)
fig.ese<- xyplot(ese~case.per.par|log.haz.lab, group=method, data=total3, xlab="# Cases per Parameter", ylab="Empirical Standard Error", type= c("l"), lty=c(1,2,1,2,1,2,1,2), col=c('black','black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0, col='grey'), main="Figure: Empirical Standard Error of Estimates", key=mykey)
fig.ese
dev.off()

pdf("ESE_3.pdf", height=9, width=9)
fig.ese.3<- xyplot(ese~case.per.par|log.haz.lab, group=method, data=subset(total3, total3$expose.rate==.3), xlab="# Cases per Parameter", ylab="Empirical Standard Error (Exposure=.3)", type= c("l"), lty=c(1,2,1,2,1,2,1,2), col=c('black', 'black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0, col='grey'), main="Figure: Empirical Standard Error of Estimates (Exposure=.3)", key=mykey)
fig.ese.3
dev.off()

pdf("ESE_5.pdf", height=9, width=9)
fig.ese.5<- xyplot(ese~case.per.par|log.haz.lab, group=method, data=subset(total3, total3$expose.rate==.5), xlab="# Cases per Parameter", ylab="Empirical Standard Error", type= c("l"), lty=c(1,2,1,2,1,2,1,2), col=c('black', 'black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0, col='grey'), main="Figure: Empirical Standard Error of Estimates (Exposure=.5)", key=mykey)
fig.ese.5
dev.off()

#########################Mean Square Error of Estimates
pdf("MSE.pdf", height=9, width=9)
fig.mse<- xyplot(mse~case.per.par|log.haz.lab, group=method, data=total3, xlab="# Cases per Parameter", ylab="Mean Square Error", type='l', lty=c(1,2,1,2,1,2,1,2), col=c('black','black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0, col='grey'), main="Figure: Mean Square Error of Estimates", key=mykey)
fig.mse
dev.off()

pdf("MSE_3.pdf", height=9, width=9)
fig.mse.3<- xyplot(mse~case.per.par|log.haz.lab, group=method, data=subset(total3, total3$expose.rate==.3), xlab="# Cases per Parameter", ylab="Mean Square Error", type='l', lty=c(1,2,1,2,1,2,1,2), col=c('black','black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0, col='grey'), main="Figure: Mean Square Error of Estimates (Exposure=.3)", key=mykey)
fig.mse.3
dev.off()

pdf("MSE_5.pdf", height=9, width=9)
fig.mse.5<- xyplot(mse~N|log.haz.lab, group=method, data=subset(total3, total3$expose.rate==.5), xlab="Sample Size", ylab="Mean Square Error", type='l', lty=c(1,2,1,2,1,2,1,2), col=c('black', 'black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0, col='grey'), main="Figure: Mean Square Error of Estimates", key=mykey)
fig.mse.5
dev.off()

#########################Median Absolute Error of Estimates
pdf("MAE.pdf", height=9, width=9)
fig.mae<- xyplot(mae~case.per.par|log.haz.lab, group=method, data=total3, xlab="# Cases per Parameter", ylab="Median Absolute Error", type='l', lty=c(1,2,1,2,1,2,1,2), col=c('black','black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0, col='grey'), main="Figure: Median Absolute Error of Estimates", key=mykey)
fig.mae
dev.off()

pdf("MAE_3.pdf", height=9, width=9)
fig.mae.3<- xyplot(mae~case.per.par|log.haz.lab, group=method, data=subset(total3, total3$expose.rate==.3), xlab="# Cases per Parameter", ylab="Median Absolute Error", type='l', lty=c(1,2,1,2,1,2,1,2), col=c('black','black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0, col='grey'), main="Figure: Median Absolute Error of Estimates (Exposure=.3)", key=mykey)
fig.mae.3
dev.off()

pdf("MAE_5.pdf", height=9, width=9)
fig.mae.5<- xyplot(mae~case.per.par|log.haz.lab, group=method, data=subset(total3, total3$expose.rate==.5), xlab="# Cases per Parameter", ylab="Median Absolute Error", type='l', lty=c(1,2,1,2,1,2,1,2), col=c('black','black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0, col='grey'), main="Figure: Median Absolute Error of Estimates (Exposure=.5)", key=mykey)
fig.mae.5
dev.off()

#########################Coverage Probability of 95% Confidence Intervals
pdf("Coverage.pdf", height=9, width=9)
fig.cover<- xyplot(Coverage~case.per.par|log.haz.lab, group=method, data=total3, xlab="# Cases per Parameter", ylab="Coverage Probability", type='l',index.cond=list(c(1,2,3,5,4)), lty=c(1,2,1,2,1,2,1,2), col=c('black','black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0.95, col='grey'), main="Figure: Coverage Probability of 95% Confidence Intervals", key=mykey)
fig.cover
dev.off()

pdf("Coverage_3.pdf", height=9, width=9)
fig.cover.3<- xyplot(Coverage~case.per.par|log.haz.lab, group=method, data=subset(total3, total3$expose.rate==.3), xlab="# Cases per Parameter", ylab="Coverage Probability", index.cond=list(c(1,2,3,5,4)), type='l', lty=c(1,2,1,2,1,2,1,2), col=c('black','black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0.95, col='grey'), main="Figure: Coverage Probability of 95% Confidence Intervals (Exposure=.3)", key=mykey)
fig.cover.3
dev.off()

pdf("Coverage_5.pdf", height=9, width=9)
fig.cover.5<- xyplot(Coverage~case.per.par|log.haz.lab, group=method, data=subset(total3, total3$expose.rate==.5), xlab="# Cases per Parameter", ylab="Coverage Probability", index.cond=list(c(1,2,3,5,4)),, type='l', lty=c(1,2,1,2,1,2,1,2), col=c('black','black', 'orange', 'orange', 'blue', 'blue', 'green', 'green'), abline=list(h=0.95, col='grey'), main="Figure: Coverage Probability of 95% Confidence Intervals (Exposure=.5)", key=mykey)
fig.cover.5
dev.off()