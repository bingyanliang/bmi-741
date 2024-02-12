rm(list=ls())

## USE THESE LINES if you don't have the packages

# install.packages("table1")
# install.packages("modeest")

library(table1)
library(modeest)

##read the data
data <- read.table("cgd.txt")
data1 <- data

head(data)

## label
data$status <- factor(data$status,
                      levels = c(0,1),
                      labels = c("censored","infection"))
data$stero <- factor(data$stero,
                     levels = c(0,1),
                     labels = c("No","Yes"))
data$proph <- factor(data$proph,
                     levels = c(0,1),
                     labels = c("No","Yes"))
data$trt <- factor(data$trt)

attach(data)
units(age) <- "years"
units(hght) <- "cm"
units(wght) <- "kg"

# original tableone
table1 <- table1(~ age + status + sex + hght + wght + inherit + stero + proph | trt, data = data )
print(table1)
table <- data.frame(table1)

## first infection

# collect data which only contains first infection and censored data
df <- data1[!duplicated(data$id), ]
# rate = number of infection/ total length of follow up
firstinfectionrate <- c("first infection rate(1/year)",
                        sum(df$status[df$trt == "placebo"])/sum(df$time[df$trt == "placebo"]),
                        sum(df$status[df$trt == "rIFN-g"])/sum(df$time[df$trt == "rIFN-g"]),
                        12*sum(df$status)/sum(df$time))
# add a row to table 1
bind <- rbind(table,firstinfectionrate)

## recurrent infection
# use original data because status should be numeric
data1 <- data.frame(data1)
# two subsets grouped by treatment method
pb_total <- subset(data1, data1$trt == "rIFN-g")
tr_total <- subset(data1, data1$trt == "placebo")

total_followup1 <- sum(pb_total$time[pb_total$status == 0])/12
total_followup2 <- sum(tr_total$time[tr_total$status == 0])/12

re_rate <- c("recurrent infection rate(1/year)",
             sum(data1$status[data1$trt == "placebo"])/total_followup1,
             sum(data1$status[data1$trt == "rIFN-g"])/total_followup2,
             12*sum(data1$status)/sum(data1$time[data1$status == 0]))
tableone <- rbind(bind, re_rate)

tableone[2,1] <- "age(years)"
tableone[11,1] <- "height(cm)"
tableone[14,1] <- "weight(kg)"

print(tableone)
View(tableone)
