rm(list = ls())

library(survival)
library(table1)

data <- read.table('invested1.txt')
df <- data

##########################################################################################
# 1. Descriptive statistics:
#    Summarize the demographic variables (gender, age, race), BMI, and medical history
#    (LVEF, MI, HF, diabetes, etc.) by treatment and overall. Be sure to also calculate the
#    (composite) event rate of death/hospitalization.
##########################################################################################
## table 1

## label
data <- na.omit(data)
data$status <- factor(data$status,
                      levels = c(0,1),
                      labels = c("censoring","death/hospitalization"))

# original tableone
table1 <- table1(~ status + gender+ age + race + dbmi + lveflt40 + priormi + priorhf + diab + renal + ischstr + pad | trtmnt, data = data )
print(table1)
table <- data.frame(table1)

# Numerator: total # of events
# (N of events is sum of status variable)
num.D <- c(sum(df$status[df$trtmnt=='HD']), 
           sum(df$status[df$trtmnt=='SD']), 
           sum(df$status))

# Denominator: total length of follow-up (year)
denom.D <- c(sum(df$time[df$trtmnt=='HD']), 
             sum(df$time[df$trtmnt=='SD']), 
             sum(df$time))/12

# Death rate
rate1 <- round(num.D/denom.D,3)
rate <- c("death/hospitalization rate(1/year)", rate1)

# add a row to table 1
tableone <- rbind(table, rate)

tableone[9,1] <- "age(years)"
#print(tableone)

###############################################################################################
# 2. Graphical analysis of subgroup effect:
#    Plot the treatment-specific (HD vs SD) Kaplan—Meier curves for hospitalization-free survival 
#    probabilities within each gender and each race. 
#    Comment on the difference in treatment effect across genders and races.
#######################################################################################################

par(mfrow=c(1,1.5))
df <- na.omit(df)
## overall
obj <- survfit(Surv(time,status)~trtmnt,data=df)
plot(obj, ylim = c(0,1), lwd=2,frame=F, lty=c(2,1),
     xlab="Time (months)",ylab="hospitalization-free survival probabilities",main="Overall",
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
legend('bottomleft', lty=2:1,c("HD","SD"),
       lwd=2,cex=1.5)

## black
obj.bl <- survfit(Surv(time,status)~trtmnt,data=df[df$race=='Black/African American',])
plot(obj.bl, xlim=c(0,80),lwd=2,frame=F, lty=c(2,1),
     xlab="Time (months)",ylab="hospitalization-free survival probabilities",main="Black/African American",
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
legend(1,0.2,lty=2:1,c("HD","SD"),
       lwd=2,cex=1.5)

## white
obj.wh<- survfit(Surv(time,status)~trtmnt,data=df[df$race=='White',])
plot(obj.wh, xlim=c(0,80),lwd=2,frame=F, lty=c(2,1),
     xlab="Time (months)",ylab="hospitalization-free survival probabilities",main="White",
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
legend(1,0.2,lty=2:1,c("HD","SD"),
       lwd=2,cex=1.5)

## Aboriginal/Native American
obj.ab<- survfit(Surv(time,status)~trtmnt,data=df[df$race=='Aboriginal/Native American',])
plot(obj.ab, xlim=c(0,80),lwd=2,frame=F, lty=c(2,1),
     xlab="Time (months)",ylab="hospitalization-free survival probabilities",main="Aboriginal/Native American",
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
legend(1,0.2,lty=2:1,c("HD","SD"),
       lwd=2,cex=1.5)

## Asian/Pacific Islander
obj.as<- survfit(Surv(time,status)~trtmnt,data=df[df$race=='Asian/Pacific Islander',])
plot(obj.as, xlim=c(0,80),lwd=2,frame=F, lty=c(2,1),
     xlab="Time (months)",ylab="hospitalization-free survival probabilities",main="Asian/Pacific Islander",
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
legend(1,0.2,lty=2:1,c("HD","SD"),
       lwd=2,cex=1.5)

## other
obj.ot<- survfit(Surv(time,status)~trtmnt,data=df[df$race=='Other',])
plot(obj.ot, xlim=c(0,80),lwd=2,frame=F, lty=c(2,1),
     xlab="Time (months)",ylab="hospitalization-free survival probabilities",main="Other",
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
legend(1,0.2,lty=2:1,c("HD","SD"),
       lwd=2,cex=1.5)

## male
obj.m <- survfit(Surv(time,status)~trtmnt,data=df[df$gender=='Male',])
plot(obj.m, xlim=c(0,80),lwd=2,frame=F, lty=c(2,1),
     xlab="Time (months)",ylab="hospitalization-free survival probabilities",main="Male",
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
legend(1,0.2,lty=2:1,c("HD","SD"),
       lwd=2,cex=1.5)

## female
obj.f<- survfit(Surv(time,status)~trtmnt,data=df[df$gender=='Female',])
plot(obj.f, xlim=c(0,80),lwd=2,frame=F, lty=c(2,1),
     xlab="Time (months)",ylab="hospitalization-free survival probabilities",main="Female",
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
legend(1,0.2,lty=2:1,c("HD","SD"),
       lwd=2,cex=1.5)


## Stratified log-rank test (by race)
survdiff(Surv(time,status)~trtmnt+ strata(race),
         data=df)
## Stratified log-rank test (by gender)
survdiff(Surv(time,status)~trtmnt+ strata(gender),
         data=df)
## Unstratified log-rank test
survdiff(Surv(time,status)~trtmnt,
         data=df)


############################################################################################
#3. Multiple regression analysis:
#   Build a Cox proportional hazards regression model to analyze the treatment effect and other risk factors for death/hospitalization. 
#   Please include the following steps:
#    a. Assess formally whether treatment effect differs significantly by gender or race (by testing interactions);
#    b. Check the proportionality assumption on the each of the covariates and take steps to address possible violations;
#    c. Check whether the linear form is appropriate for age and BMI; make proper transformation/categorization in case of non-linearity;
#    d. Interpret the analysis results of your final model, e.g., what risk factors are most strongly or significantly associated with 
#       the primary endpoint and to what degree.
#############################################################################################

# a. Assess formally whether treatment effect differs significantly by gender or race (by testing interactions);

#fit the Cox proportional hazards model
obj <- coxph(Surv(time,status)~age+factor(gender)+factor(race)+dbmi+factor(lveflt40)+factor(priormi)
             +factor(priorhf)+factor(diab)+factor(renal)+factor(ischstr)+factor(pad),data=df)
summary(obj)


# b. Check the proportionality assumption on the each of the covariates and take steps to address possible violations;

#Schoelfeld residuals
#produce proportionality test results
sch <- cox.zph(obj) 
print(sch) 
#Plot scaled Schoelfeld residuals for each covariate
par(par(mar = c(4, 4, 2, 2)), mfrow=c(3,3))
plot(sch,xlab="Time (months)",lwd=2,cex.lab=1.2,cex.axis=1.2)

obj.stra <- coxph(Surv(time,status)~age+factor(race)+dbmi+factor(lveflt40)+factor(priormi)+factor(priorhf)+factor(diab)
                  +factor(renal)+factor(ischstr)+factor(pad)+strata(gender),data=df)
sch.stra <- cox.zph(obj.stra)
print(sch.stra)

#fit the Cox proportional hazards model
obj.stra <- coxph(Surv(time,status)~age+factor(race)+dbmi+factor(lveflt40)+factor(priorhf)+factor(diab)
                  +factor(renal)+factor(ischstr)+factor(pad)+strata(gender)+strata(priormi),data=df)
sch.stra <- cox.zph(obj.stra)
print(sch.stra)

# c. Check whether the linear form is appropriate for age and BMI; 
# make proper transformation/categorization in case of non-linearity;

## Martingale residuals
mart_resid <- resid(obj,type='martingale')

#plot the martingale residuals against age and bmi

## Age
par(mfrow=c(2,2))
plot(df$age, mart_resid,
     xlab="Age (years)", ylab="Martingale residuals",
     main='Age',cex.lab=1.5)
lines(lowess(df$age, mart_resid),lwd=2)
abline(0,0,lty=3,lwd=2)

## bmi
plot(log(df$dbmi), mart_resid,
     xlab="Tumor size (mm)", ylab="Martingale residuals",
     main='Tumor size',cex.lab=1.2)
lines(lowess(log(df$dbmi), mart_resid),lwd=2)
abline(0,0,lty=3,lwd=2)

# To address non-linear age
# categorize age in agec
# age<=40: agec=1
# 40<age<=60: agec=2
# 60<age<=80: agec=3
# age>80: agec=4
df$agec <- (df$age<=40)+2*(df$age>40&df$age<=60)+3*(df$age>60&df$age<=80)+4*(df$age>80)
df$bmic <- (df$dbmi<=30)+2*(df$dbmi>30&df$dbmi<=50)+3*(df$dbmi>50&df$dbmi<=70)+4*(df$dbmi>70)

#re-fit the model with agec and bmic
obj.stra.final <- coxph(Surv(time,status)~factor(agec)+factor(race)+factor(bmic)+factor(lveflt40)+factor(priorhf)+factor(diab)+factor(renal)+factor(ischstr)+factor(pad)+strata(gender)+strata(priormi),data=df)


#plot the estimated HRs for agec (Figure 4.6)
# and confidence intervals
final.sum <- summary(obj.stra.final)
final.sum


# Plot the age-group-specific HR and confidence
# intervals from the re-fitted model
ci.table=final.sum$conf.int
hr=ci.table[3:5,1]
hr.low=ci.table[3:5,3]
hr.up=ci.table[3:5,4]

par(mfrow=c(1,1))
plot(1:4,c(1,hr),ylim=c(0,1.5),frame=F,xaxt='n',
     xlab="Age (years)", ylab="Hazard ratio",pch=19,cex=1.5,cex.lab=1.2,
     cex.axis=1.2)
axis(1, at=c(1,2,3,4),labels=c("(20, 40]","(40, 60]","(60, 80]","(80, 100]"),cex.axis=1.5)
# horizontal error bars
arrows(2:4, hr.low, 2:4, hr.up, length=0.05, angle=90, code=3, lwd=2)
lines(1:4,c(1,hr),lty=3,lwd=2)

# Plot the bmi-group-specific HR and confidence
# intervals from the re-fitted model

hr=ci.table[10:12,1]
hr.low=ci.table[10:12,3]
hr.up=ci.table[10:12,4]

par(mfrow=c(1,1))
plot(1:4,c(1,hr),ylim=c(0,3),frame=F,xaxt='n',
     xlab="bmi", ylab="Hazard ratio",pch=19,cex=1.5,cex.lab=1.2,
     cex.axis=1.2)
axis(1, at=c(1,2,3,4),labels=c("(10, 30]","(30, 50]","(50, 70]","(70, 90]"),cex.axis=1.2)
# horizontal error bars
arrows(2:4, hr.low, 2:4, hr.up, length=0.05, angle=90, code=3,lwd=2)
lines(1:4,c(1,hr),lty=3,lwd=2)

# d. Interpret the analysis results of your final model, e.g., what risk factors are most strongly or significantly associated with 
#   the primary endpoint and to what degree.




