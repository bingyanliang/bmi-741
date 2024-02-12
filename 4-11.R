rm(list=ls())
library(survival)

gbc <- read.table("gbc.txt")

o <- order(gbc$id,gbc$time)
gbc <- gbc[o,]
#get the first row for each id
data.CE <- gbc[!duplicated(gbc$id),]

#set status=1 if status==2 or 1
data.CE$status <- (data.CE$status>0)+0

#fit the Cox proportional hazards model
obj <- coxph(Surv(time,status)~factor(hormone)+factor(meno)+age+size+factor(grade)
             +prog+estrg,data=data.CE)
#summarize results
summary(obj)

### Get the Breslow estimator for baseline
### cumulative hazard function
Lambda0 <- basehaz(obj,centered=F)
coef <- obj$coefficients
z1 <- c(0, 0, median(data.CE$age), median(data.CE$size),1,0,median(data.CE$prog),median(data.CE$estrg))
z2 <- c(1, 0, median(data.CE$age), median(data.CE$size),1,0,median(data.CE$prog),median(data.CE$estrg))
s1 <- exp(-exp(sum(coef*z1))*Lambda0$hazard)
s2 <- exp(-exp(sum(coef*z2))*Lambda0$hazard)
plot(Lambda0$time, s1, type = 'l', col=1, xlab = 'time', ylab = 'survival')
lines(Lambda0$time, s2, type = 'l', col =2)
