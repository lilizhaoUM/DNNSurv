rm(list=ls())
library(keras)
install_keras()
library(pseudo)
library(survivalROC)
library(survival)
library(survcomp)
library(survAUC)

# simulate survival data from a cox model with one covariate 
set.seed(2412)
pickTime=c(0.7,1.7,3.2,5.3,8.3)
n=2000
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
x_train <- scale(x_train, center = mean, scale = std)

# data normalization
x_test=x[(n/2+1):n]
surv_test <- time[(n/2+1):n]
cen_test =delta[(n/2+1):n]
x_test<- scale(x_test, center = mean, scale = std)


pseudoCond  = getPseudoConditional(surv_train, cen_train, pickTime)

x <- x_train[pseudoCond$id,]

# create dummy variables for the time points
smatrix=model.matrix(~as.factor(pseudoCond$s)+0)
x_train.all <- cbind(x, smatrix)

y_train.all <- pseudoCond$pseudost

model = pseudoDNN.train(x_train.all, y_train.all)


# format the test data 
x_test.all=do.call(rbind, replicate(length(pickTime), x_test, simplify=FALSE))

s_test=rep(s,each=nrow(x_test))
smatrix.test=model.matrix(~as.factor(s_test)+0)
x_test.all=cbind(x_test.all,smatrix.test)

# predict test data
ypred.con <- pseudoDNN.predict(model, x_test.all)

# format the conditional survival probability at  each time point
ypred.con <- matrix(ypred.all, nrow=nrow(x_test))

# obtain the marginal survival probability by multiple series of conditional probabilities
ypred <- lapply(1:length(s), function(i) apply(ypred.con[,1:i, drop=FALSE], 1, prod))
ypredsur <- Reduce(cbind, ypred)


