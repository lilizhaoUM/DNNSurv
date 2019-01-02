rm(list=ls())
library(keras)
install_keras()
library(pseudo)
library(survivalROC)
library(survival)
library(survcomp)
library(survAUC)

# simulate survival data from a cox model with one covariate 
# the censoring distribution depends on the covariate

set.seed(2412)
pickTime=c(0.7,1.7,3.2,5.3,8.3)
x=rnorm(n,0,1)
c0=0.1
times=rexp(n,c0*exp(1*x))  
time.censor=rexp(n,c0*exp(1*x)) 
summary(time.censor)
delta=ifelse(times<time.censor, 1, 0)
time <- ifelse(times<time.censor, times, time.censor)

x_train  <- x[1:(n/2)]
surv_train <- time[1:(n/2)]
cen_train <- delta[1:(n/2)]

# data normalization
mean <- apply(as.matrix(x_train), 2, mean)
std <- apply(as.matrix(x_train), 2, sd)
xs_train <- scale(x_train, center = mean, scale = std)

# data normalization
x_test=x[(n/2+1):n]
surv_test <- time[(n/2+1):n]
cen_test =delta[(n/2+1):n]
xs_test<- scale(x_test, center = mean, scale = std)

# obtain the risk probability for each pickTime
ypredw <- predictSurvprob(surv_train, cen_train, weight=1, pickTime, xs_train, xs_test) 

# obtain the c-index (first column) and Brier score (last column)
evaluateCondPseudo(1-ypredw, surv_train, cen_train, pickTime, surv_test, cen_test)


# IPCW method
weight=compute_weight(surv_train,cen_train,x_train)
ypredw <- predictSurvprob(surv_train, cen_train, weight=weight, pickTime, xs_train, xs_test) 


