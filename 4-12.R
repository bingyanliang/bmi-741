rm(list = ls())
library(survival)

###################################################
# Section 4.2.4 (Figure 4.2)                
# Cox proportional hazards model on GBC study data
###################################################
brca <- read.table("BrCa.txt")

#subset to first event data
#Sort the data by time within each id
o <- order(brca$id,brca$time)
brca <- brca[o,]

#get the first row for each id
brca.CE <- brca[!duplicated(brca$id),]

#set status=1 if status==2 or 1
brca.CE$status <- (brca.CE$status>0)+0


#fit the Cox proportional hazards model
obj <- coxph(Surv(time,status)~age+factor(meno)+factor(size)+factor(grade)+nodes+prog+estrg
             +factor(horm)+factor(chem),data=brca.CE)
#summarize results
round(summary(obj)$coef, 3)


############################
# Residual analysis
# (Figures 4.3--4.6)
############################


#####################################

## First get the Cox-Snell residuals.;
## The default residuals of coxph in R are the martingale residuals.
## resid(obj,type=c("martingale", "deviance", "score", "schoenfeld",
##                   "dfbeta", "dfbetas", "scaledsch","partial"))



## Use relationship between cox-snell and martingal
## residuals
coxsnellres <- brca.CE$status-resid(obj,type="martingale")
## Then use N-A method to estimate the cumulative 
## hazard function for residuals;
fit <- survfit(Surv(coxsnellres,brca.CE$status)~1)
Htilde <- cumsum(fit$n.event/fit$n.risk)
plot(log(fit$time),log(Htilde),cex.axis=1.5,cex.lab=1.5,lwd=1.5,
     main="Cox-Snell residual plot",xlab="log t", cex.main=1.5,
     ylab="log cumulative hazard")
abline(0,1,lty=2,lwd=1.5)

#Schoelfeld residuals
#produce proportionality test results
sch <- cox.zph(obj) 
print(sch) 
#Plot scaled Schoelfeld residuals for each covariate
par(par(mar = c(4, 4, 2, 2)), mfrow=c(3,3))
plot(sch,xlab="Time (years)",lwd=2,cex.lab=1.2,cex.axis=1.2)



# To address non-proportionality of tumor grade
# re-fit the Cox proportional hazards model
# stratified by tumor grade
obj.stra <- coxph(Surv(time,status)~age+factor(meno)+factor(size)+strata(grade)+nodes+prog+estrg
                  +factor(horm)+factor(chem),data=brca.CE)
#Schoelfeld
#produce proportionality test results
sch.stra <- cox.zph(obj.stra) 
print(sch.stra) 
par(mfrow=c(5,2))
plot(sch.stra,xlab="Time (years)",lwd=2,cex.lab=1.5)

## Martingale residuals
mart_resid <- resid(obj,type='martingale')

#plot the martingale residuals against
# the fours quantitive covariates:
# age, nodes, progesterone and estrogen receptor levels

## age
par(mfrow=c(2,2))
plot(brca.CE$age, mart_resid,
     xlab="age (years)", ylab="Martingale residuals",
     main='Age',cex.lab=1.2,cex.axis=1.2)
lines(lowess(brca.CE$age, mart_resid),lwd=2)
abline(0,0,lty=3,lwd=2)

## nodes
plot(brca.CE$nodes, mart_resid,
     xlab="Number of positive lymph nodes", ylab="Martingale residuals",
     main='Number of positive lymph nodes',cex.lab=1.2,cex.axis=1.2)
lines(lowess(brca.CE$size, mart_resid),lwd=2)
abline(0,0,lty=3,lwd=2)

## Progesterone
plot(brca.CE$prog, mart_resid,
     xlab="Progesterone receptor (fmol/mg)", ylab="Martingale Residuals",
     main='Progesterone',cex.lab=1.2,cex.axis=1.2)
lines(lowess(brca.CE$prog, mart_resid),lwd=2)
abline(0,0,lty=3,lwd=2)


## Estrogen
plot(brca.CE$estrg, mart_resid,
     xlab="Estrogen receptor (fmol/mg)", ylab="Martingale Residuals",
     main='Estrogen',cex.lab=1.2,cex.axis=1.2)
lines(lowess(brca.CE$estrg, mart_resid),lwd=2)
abline(0,0,lty=3,lwd=2)



