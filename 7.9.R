rm(list = ls())

library(survival)
library(IntCens)


## BMA HIV study
## Bangkok Metropolitan Administration HIV Study
data <- read.table("bam.txt")


# get the response data for ICSurvICM()
delta <- data$delta
gamma <- data$gamma
n <- nrow(data)
U <- data$U
V <- data$V


# Fit a proportional odds model
obj.PO <- ICSurvICM(delta,gamma,U,V,Z=data[,5:9],model="PO")
print.ICSurv(obj.PO)

# Table 7.1
# construct a table for odd ratio and
# 95% confidence intervals

#regression parameter
beta <- obj.PO$beta
se <- sqrt(diag(obj.PO$var))
c1 <- round(exp(beta),2)
c2 <- paste0("(",round(exp(beta-1.96*se),2),", ",
             round(exp(beta+1.96*se),2),")")
noquote(cbind(c1,c2))


