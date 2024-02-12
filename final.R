rm(list = ls())

library(survival)
library(tidyverse)

# Read in data
invested <- read.table("invested2.txt", header = TRUE)
invested <- na.omit(invested)
invested <- invested[invested$time>0, ]

## 1.1
# Create table of event rates by treatment group and overall for death
data <- invested$time[invested$status==0 | invested$status==2]
overall_time <- sum(data)
invested2 <- invested[!duplicated(invested$patid),]
sum1 <- sum(invested$status==2)
sum2 <- sum(invested2$status==2)
sum3 <- sum1 - sum2
death_table <- invested %>%
  filter(status == 2) %>%
  group_by(trtmnt) %>%
  summarize(overall = sum(status == 2)/overall_time,
            before_hosp = sum2/overall_time,
            after_hosp = sum3/overall_time,
            n = n())

# Create table of event rates by treatment group and overall for first hospitalization
hosp_table <- invested2 %>%
  filter(status == 1) %>%
  group_by(trtmnt) %>%
  summarize(overall = sum(status == 1)/overall_time,
            n = n())

# Create table of event rates by treatment group and overall for recurrent hospitalizations
library(dplyr)

rec_hosp <- invested %>%
  group_by(trtmnt) %>%
  summarise(event_rate = sum(status == 1)/overall_time) %>%
  bind_rows(invested %>%
              summarise(trtmnt = "Overall",
                        event_rate = sum(status == 1)/overall_time))



# Print tables
death_table
hosp_table
rec_hosp


##1.2
# Calculate number of deaths/hospitalizations lost by only looking at time to first event
lost <- length(invested$status!=0) - length(invested2$status!=0)
lost

## 1.3
# Create histogram of number of hospitalizations by treatment group
ggplot(invested, aes(x = factor(patid), fill = factor(status == 1))) +
  facet_wrap(~trtmnt, ncol = 2) +
  geom_bar() +
  labs(x = "Patient", y = "Number of Hospitalizations") +
  scale_fill_manual(values = c("#0072B2", "#F0E442")) +
  theme_minimal()


## 2.1
library(cmprsk)

df <- invested

# recode trtmnt as factor
df$trtmnt <- factor(df$trtmnt, levels = c("SD", "HD"))

# Convert character variables to factors
df$trtmnt <- as.numeric(df$trtmnt)
df$gender <- as.numeric(df$gender)
df$race <- as.numeric(df$race)
df$lveflt40 <- as.numeric(df$lveflt40)
df$priormi <- as.numeric(df$priormi)
df$priorhf <- as.numeric(df$priorhf)
df$diab <- as.numeric(df$diab)
df$renal <- as.numeric(df$renal)
df$ischstr <- as.numeric(df$ischstr)
df$pad <- as.numeric(df$pad)

# Fit the competing risks model
obj.fg <- crr(df$time, df$status, df[,5:15], failcode=1)

# Estimate cumulative incidence
cuminc.obj <- cuminc(df$status, df$time, group = df$gender, rho=0)

# Plot cumulative incidence
plot(cuminc.obj, xlab = "Time (months)", ylab = "Cumulative incidence")

# Test for treatment effect
print(summary(obj.fg))


## 2.2.1
# Fit a Cox proportional hazards model with shared frailty for recurrent hospitalizations
frailty_model <- coxph(Surv(time, status == 1) ~ trtmnt + gender + age + race + dbmi + 
                         lveflt40 + priormi + priorhf + diab + renal + ischstr + pad + frailty(patid), 
                       data = invested)

# Print the summary of the model
summary(frailty_model)


## 2.2.2
# Fit a Cox proportional hazards model for overall survival
survival_model <- coxph(Surv(time, status != 0) ~ trtmnt + gender + age + race + dbmi + 
                          lveflt40 + priormi + priorhf + diab + renal + ischstr + pad, 
                        data = invested)

# Print the summary of the model
summary(survival_model)


## 2.2.3
# The direction and magnitude of the covariate effects can be interpreted from the coefficient 
# estimates and corresponding p-values in the summary output of the models. 
# Positive coefficients indicate higher risk, while negative coefficients indicate lower risk. 
# A p-value less than 0.05 indicates statistical significance.


## 2.3
library(rmt)
library(dplyr)


# Create a composite endpoint variable
invested <- invested %>%
  mutate(comp = if_else(status == 2, 1, 2))

# Specify the restricting time tau and the maximum number K of recurrent events to consider
tau <- 12 # months
K <- 2

# Compute the RMTIF of treatment using the composite endpoint with death prioritized over recurrent hospitalizations
rmt <- rmtfit(invested$patid, invested$time, invested$status, invested$trtmnt, type = "multistate")
summary(rmt)


#3.1
data <- read.table('invested2.txt')
library("glmnet")
## Loading required package: Matrix
## Loaded glmnet 4.1-7
library(rpart)
library(rpart.plot)
o <- order(data$patid,data$time)
dat <- data[o,]
data.CE <- dat[!duplicated(dat$patid),]
data.CE <- data.CE[data.CE$time>0,]
data.CE$status <- (data.CE$status>0)+0
set.seed(1234)
n <- nrow(data.CE)
ind <- sample(1:n)[1:4000]
train <- data.CE[ind,]
test <- data.CE[-ind,]
set.seed(12345)
obj <- rpart(Surv(time, status) ~ factor(trtmnt)+gender+age+factor(race)+dbmi+factor(lveflt40)+factor(priormi) + factor(priorhf) + 
               factor(diab) + factor(renal) + factor(ischstr) + factor(pad), 
             data = data.CE )
cptable <- obj$cptable
# complexity parameter values
CP <- cptable[,1]
             # obtain the optimal parameter
cp.opt <- CP[which.min(cptable[,4])]
             # Prune the tree
fit <- prune(obj, cp = cp.opt)
rpart.plot(fit)


## 3.2
newpat1 <- data[1,]
newpat1$trtmnt <- 'SD'
newpat1$gender <- 'Female'
newpat1$age <- 60
newpat1$race <- 'White'
newpat1$dbmi <- 25
newpat1$lveflt40 <- 'Yes'
newpat1$priormi <- 'Yes'
newpat1$priorhf <- 'Yes'
newpat1$diab <- 'Yes'
newpat1$renal <- 'Yes'
newpat1$ischstr <- 'Yes'
newpat1$pad <- 'Yes'
newpat2 <- newpat1
newpat2$trtmnt <- 'HD'
# Extract predicted terminal node for each observation
train$node <- predict(fit, newdata = train, type = "vector")

# Get the KM estimates for the outcome in each terminal node
km <- survfit(Surv(time, status) ~ node, data = train)
tmp <- summary(km)
tmp.strata <- as.integer(sub(".*=", "", tmp$strata))
tmp.t <- tmp$time
tmp.surv <- tmp$surv
# Number of terminal nodes
TN <- unique(tmp.strata)
N <- length(TN)
# Combine the predicted survival rates together,
# as functions of t
t <- sort(unique(tmp.t))
m <- length(t)
fitted_surv=matrix(NA,m,N)
for (j in 1:m){
  tj <- t[j]
  for (k in 1:N){
    tk <- c(0,tmp.t[tmp.strata==TN[k]])
    survk <- c(1,tmp.surv[tmp.strata==TN[k]])
    fitted_surv[j,k] <- survk[sum(tk<=tj)]
  }
}
# Get the terminal node prediction
# for the test data
library(treeClust)
test_term1 <- rpart.predict.leaves(fit, newpat1)
n <- length(test_term1)
St_tree1 <- matrix(NA,n,m)
for (k in 1:N){
  ind <- which(test_term1==TN[k])
  St_tree1[ind,] <- matrix(fitted_surv[,k], nrow=length(ind),
                           ncol=m, byrow=TRUE)
}
plot(t, St_tree1)
