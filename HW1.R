


## Read in the CGD  data

data <- read.table("Chronic Granulomatous Disease Study//cgd.txt")
head(data)


###########################################################
# Table 1 Patient characteristics for the CGD study
###########################################################


##A function calculating median (IQR) by 
## binary group
## Input: y=quantitative variable
##        trt=binary group variable
##        decp=number of decimal points
## Output: a row vector containing median (IQR)
##    by the two levels of trt and overall 

Mean.IQR.by.trt=function(y,trt,decp=1){
  groups=sort(unique(trt))
  all=quantile(y)
  g1=quantile(y[trt==groups[1]])
  g2=quantile(y[trt==groups[2]])
  
  result=matrix(NA,1,3)
  colnames(result)=c(groups,"Overall")
  result[1,1]=paste0(round(g1[3],decp)," (",round(g1[2],decp),", ",round(g1[4],decp),")")
  result[1,2]=paste0(round(g2[3],decp)," (",round(g2[2],decp),", ",round(g2[4],decp),")")
  result[1,3]=paste0(round(all[3],decp)," (",round(all[2],decp),", ",round(all[4],decp),")")
  return(result)
}


##A function calculating N (%) by 
## binary group
## Input: x=categorical variable with p levels
##        trt=binary group variable
##        decp=number of decimal points of %
## Output: a px3 matrix containing N (%) for each level of x
##    by the two levels of trt and overall 

N.prct.by.trt=function(x,trt,decp=1){
  groups=sort(unique(trt))
  x.levels=sort(unique(x))
  p=length(x.levels)
  n=length(x)
  n1=length(x[trt==groups[1]])
  n2=length(x[trt==groups[2]])
  
  result=matrix(NA,p,3)
  colnames(result)=c(groups,"Overall")
  rownames(result)=x.levels
  
  for (i in 1:p){
    n1i=sum(x[trt==groups[1]]==x.levels[i])
    n2i=sum(x[trt==groups[2]]==x.levels[i])
    ni=sum(x==x.levels[i])
    
    
  result[i,1]=paste0(n1i," (",round(n1i/n1*100,decp),"%)")
  result[i,2]=paste0(n2i," (",round(n2i/n2*100,decp),"%)")
  result[i,3]=paste0(ni," (",round(ni/n*100,decp),"%)")
  }
  
  
  return(result)
}


##Baseline characteristics by treatment arm:
# sex
# Sex of each patient(male, female)
# 
# age
# Age of each patient at study entry, in years
# 
# height
# Height of each patient at study entry, in cm
# 
# weight
# Weight of each patient at study entry, in kg
# 
# inherit
# Pattern of inheritance (autosomal recessive, X-linked)
# 
# steroids
# Using corticosteroids at times of study centry(1=Yes, 0=No)
# 
# propylac
# Using prophylactic antibiotics at time of study entry(1=Yes, 0=No)


####################################################
# To calculate the baseline characteristics,
# first subset to one record per patient.
# We do this by getting the first record 
#####################################################

o <- order(data$id,data$time)
dat <- data[!duplicated(data$id),]
n <- nrow(dat)
#n=128: number of subjects

table1 <- rbind(
  N.prct.by.trt(x=dat$sex,trt=dat$trt),
  Mean.IQR.by.trt(y=dat$age,trt=dat$trt),
  Mean.IQR.by.trt(y=dat$hght,trt=dat$trt),
  Mean.IQR.by.trt(y=dat$wght,trt=dat$trt),
  N.prct.by.trt(x=dat$inherit,trt=dat$trt),
  N.prct.by.trt(x=dat$stero,trt=dat$trt),
  N.prct.by.trt(x=dat$proph,trt=dat$trt)
)


noquote(table1)

#############################################
##Calculate the *first event* rates (per year)
# use the dataset "dat", since it contains
# the first record for each patient
############################################


#Numerator: total # of events
# note that event is coded as status=2
num.FE <- c(sum(dat$status[dat$trt=="placebo"]>0), 
        sum(dat$status[dat$trt=="rIFN-g"]>0), 
        sum(dat$status>0))

#Demoninator: total length of follow-up (year)
denom.FE <- c(sum(dat$time[dat$trt=="placebo"]), 
        sum(dat$time[dat$trt=="rIFN-g"]), 
        sum(dat$time))/12

#death rate
round(num.FE/denom.FE,3)

#############################
#Recurrent event rate       #
#############################

# the numerator is easy: simply the numbers of records
# with non-zero status
num.rec <- c(sum(data$status[data$trt=="placebo"]>0), 
         sum(data$status[data$trt=="rIFN-g"]>0), 
         sum(data$status>0))

# Denominator should be the sum of the longest
# follow-up for each patient

#sort the data by descending order of time
o <- order(data$id,data$time,decreasing = TRUE)
data <- data[o,]
# get the first record, which contains
# the overall length of follow-up
# for each patient
dat.last <- data[!duplicated(data$id),]

denom.rec <- c(sum(dat.last$time[dat.last$trt=="placebo"]), 
            sum(dat.last$time[dat.last$trt=="rIFN-g"]), 
            sum(dat.last$time))/12

#recurrent event rate
round(num.rec/denom.rec,3)




