rm(list=ls())
# Part One: the baseline(table1)

# read data
data <- read.table("smoking.txt")

# install.packages("table1")
# install.packages("modeest")

library(table1)
library(modeest)

## label
data$status <- factor(data$status,
                      levels = c(0,1),
                      labels = c("censored","return"))

# original tableone
table1 <- table1(~ age + sex + race + empl + years + level + pa + nosmk | trt, data = data )

# returnrate = number of return/ total length of follow up
smoking <- read.table("smoking.txt")
table <- data.frame(table1)

smoking$trt <- factor(smoking$trt)
c1 <- sum(smoking$status[smoking$trt == "combination"])/sum(smoking$time[smoking$trt == "combination"])
c2 <- sum(smoking$status[smoking$trt == "patchOnly"])/sum(smoking$time[smoking$trt == "patchOnly"])
c3 <- sum(smoking$status)/sum(smoking$time)
returnrate <- c("return rate(1/year)",365.25*c1,365.25*c2,365.25*c3)

# add a row to table 1
tableone <- rbind(table,returnrate)
tableone[2,1] <- "age(years)"
tableone[13,1] <- "employment status"
tableone[14,1] <- "  full time"
tableone[16,1] <- "  part time"
tableone[17,1] <- "smoking years"
tableone[23,1] <- "past attempts"
tableone[26,1] <- "nonsmoke(years)"

print(tableone)

# Part Two: Analyzing
################################################
# the difference is number of black people in two groups
# Revise Professor's sample code of Figure 3.5  
################################################
library(survival)
par(mfrow=c(1,1.5))

## overall
obj <- survfit(Surv(time,status)~trt,data=smoking)
plot(obj, ylim = c(0,1), lwd=2,frame=F, lty=c(2,1),
     xlab="Time (days)",ylab="Relaspe-free survival probabilities",main="Overall",
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
legend('bottomleft', lty=2:1,c("combination","patchOnly"),
       lwd=2,cex=1.5)

## Stratified log-rank test (by race)
survdiff(Surv(time,status)~trt+ strata(race),
         data=smoking)
## Unstratified log-rank test
survdiff(Surv(time,status)~trt,
         data=smoking)


