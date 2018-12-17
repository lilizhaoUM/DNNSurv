rm(list=ls())
library(keras)
install_keras()
library(pseudo)
library(survivalROC)
library(survival)
library(survcomp)
library(survAUC)

#----------------------EvaluateCondPseudo-------------------------------------
#-ypred is risk probabilities from the predictio model
#-pickTime is the time points used to divivide the time intervals
#----------------------------------------------------------------------------
evaluateCondPseudo <- function(ypred, surv_train, cen_train, pickTime,
                               surv_test, cen_test){
  
  measures.all <- matrix(0, length(pickTime), 4)
  for(i in 1:length(pickTime)){
    measures.all[i,] <- unlist(evaluate(surv_train, cen_train, ypred[,i], 
                                        surv_test, cen_test, pickTime[i]))
  }
  measures.all
}
#----------------------------------------------------------------------------
# evaluate the performance of survival prediction using different metrics
# ypred is a risk probability
#----------------------------------------------------------------------------
evaluate <-function(surv_train, cen_train, ypred,
                    surv_test, cen_test, pt){
  c_index <- concordance.index(x=ypred, surv.time=surv_test,
                               surv.event=cen_test, method="noether")$c.index
  
  auc <- survivalROC.C(Stime = surv_test,
                       status = cen_test,
                       marker =  ypred,
                       predict.time = pt)$AUC
  
  surv.rsp <- Surv(surv_train, cen_train)
  surv.rsp.new <- Surv(surv_test, cen_test)
  unoc <- UnoC(surv.rsp, surv.rsp.new, ypred)
  
  brier <- BrierScore(surv_test, cen_test, pt, 1-ypred, type = "kaplan")
  
  list(c_index=c_index, unoc=unoc, auc=auc, brier=brier)

}

#----------------------------------------------------------------------------
# calcute censoring distribution
#----------------------------------------------------------------------------
Ghat.FUN <- function(time, status,  new.time, type = "kaplan"){
  ranked.new.time <- rank(new.time)
  summary(survfit(Surv(time, status) ~ 1, se.fit = FALSE, type = type), sort(new.time))$surv[ranked.new.time]
}

#----------------------------------------------------------------------------
# calcuate the Brier score using time, status, survprob from the prediction model
#----------------------------------------------------------------------------
BrierScore <- function(time, status, tau, survprob, type = "kaplan"){
  
  # remove sub who is censored before tau
  rm <- (1-status)*(time < tau) 
  sub.time=time[!rm]  
  sub.status=status[!rm]
  sub.survprob=survprob[!rm]
  numerator   <- (sub.survprob)^2*sub.status * (sub.time <= tau) + (1-sub.survprob)^2*(sub.time >= tau)
  min.time    <- pmin(sub.time, tau)
  denominator <- Ghat.FUN(time, 1 - status, min.time, type = type)
  result      <- numerator/denominator
  result[is.na(result)] <- 0
  return(mean(result))
}


#----------------------------------------------------------------------------
#-two-layer deep neural network model in keras
#----------------------------------------------------------------------------
pseudoDNN.train <- function(x_train, y_train){
  
  model <- keras_model_sequential() %>%
    #use tanh instead of relu for activation function
    layer_dense(units = 8, kernel_regularizer = regularizer_l2(0.0001), activation = "tanh",
                input_shape = dim(x_train)[[2]]) %>%
    layer_dense(units = 4, kernel_regularizer = regularizer_l2(0.01),
                activation = "tanh") %>%
    layer_dense(units = 1, activation='sigmoid')
  
  #use adam instead of rmsprop for optimizer
  model %>% compile(
    optimizer = optimizer_adam(lr = 0.0025),
    loss = "mse",
    metrics = c("mae")
  )
  model %>% fit(x_train, y_train,
                epochs = 30, batch_size = 64,
                verbose = 0)
  
  model
}

#----------------------------------------------------------------------------
#prediction based on the above two-layer keras model
#----------------------------------------------------------------------------
pseudoDNN.predict <- function(model, x_test){
  ypred <- model %>% predict(x_test)
  ypred
}

#----------------------------------------------------------------------------
# -compute weights for the ipcw method (which is 1/surv censoring probability)
# -output is n*ns matrix: ns is the number of unique time points 
#----------------------------------------------------------------------------
compute_weight<- function(time,delta,cov){
  ## preparing the data
  n=length(time)
  dat <- data.frame(id=1:n,t=time,delta=delta,cov)
  
  # sort in time, if tied, put events before censoring!!
  dat <- dat[order(dat$t,-dat$delta),]
  cov <- dat[,-c(1,2,3)]
  t=dat$t
  delta=dat$delta
  # coef for censoring
  censor_model=coxph(Surv(t,1-delta)~., data=dat[,-1])
  beta_censor=coef(censor_model)
  
  n=length(t)
  s=sort(unique(t))
  
  ns=length(s)  # the number of intervals
  C <- do.call(cbind, lapply(1:ns, function(j)  (s[j] == t)*(delta == 0)))
  Y <- do.call(cbind, lapply(1:ns, function(j) ifelse(s[j] <= t, 1, 0)))
  
  lp=as.vector(as.matrix(cov,nrow=n)%*%beta_censor)
  denominator_weight=colSums(Y*exp(lp))
  numerator_weight=outer(colSums(C),exp(lp),  "*") 
  Haz=apply(numerator_weight/denominator_weight,2,cumsum)
  survprob=t(exp(-Haz))  # n*ns dimension
  trunc=max(1e-20,min(survprob[survprob>0])) # avoid close to zero prob
  weight=ifelse(survprob < trunc, 1/trunc, 1/survprob)
  weight <- weight[order(dat$id),]                             #back to original order
  return(weight)
}



##############################################################################################
# get pseudo conditional survival probabilities 
# weight allows IPCW calcuations. If weight=1, it is the non-IPCW method
# output is a n*length(tau)  matrix
#######################################################################################
getPseudoConditional<-function(time,delta,weight,tau){
  
  ## preparing the data
  n=length(time)
  ntau=length(tau)
  pseudo <- data.frame(id=1:n,t=time,delta=delta,weight)
  
  # sort in time, if tied, put events before censoring
  pseudo <- pseudo[order(pseudo$t,-pseudo$delta),]
  
  t=pseudo$t
  delta=pseudo$delta
  weight=pseudo[,-c(1,2,3)]
  
  s=sort(unique(t))
  ns=length(s)  # the number of intervals
  D <- do.call(cbind, lapply(1:ns, function(j)  (s[j] == t)*(delta == 1)))
  Y <- do.call(cbind, lapply(1:ns, function(j) ifelse(s[j] <= t, 1, 0)))
  
  Yw=Y*weight
  Dw=D*weight
  denominator=colSums(Yw)
  numerator=colSums(Dw)
  IPCW_CH=cumsum(numerator/denominator)
  IPCW_CH.teval=approx(s, IPCW_CH, xout = tau,method="constant")$y
  IPCW_surv.teval=exp(-IPCW_CH.teval)
  
  Denominator=matrix(denominator,n,ns,byrow=TRUE)-Yw
  Numerator=matrix(numerator,n,ns,byrow=TRUE)-Dw
  IPCWi_CH=apply(Numerator/Denominator,1,cumsum)
  IPCWi_CH.teval=t(sapply(1:n,function(j) approx(s, IPCWi_CH[,j], xout = tau,method="constant")$y))
  
  IPCW_survi.teval=matrix(exp(-IPCWi_CH.teval),nrow=n)
  mat.tmpi=cbind(rep(1,n),IPCW_survi.teval[,-ntau])
  con.IPCW_survi.teval=IPCW_survi.teval/mat.tmpi
  
  surv.all.mat=matrix(IPCW_surv.teval,nrow=n,ncol=ntau,byrow=T)
  mat.tmp=cbind(rep(1,n),surv.all.mat[,-ntau])
  con.surv.all.mat=surv.all.mat/mat.tmp
  
  # adjust n at risk 
  s.tau=c(0,tau[-ntau])
  R <- do.call(cbind, lapply(1:length(s.tau), function(j) ifelse(t>s.tau[j], 1, 0)))
  nr=colSums(R)
  nr.mat=matrix(nr,nrow=n,ncol=ntau,byrow=TRUE)
  con_pseudo_survprob=nr.mat*con.surv.all.mat-(nr.mat-1)*con.IPCW_survi.teval
  
  # back to original order
  out <- NULL
  out$tau <- tau
  out$conditional_pseudost <- as.matrix(con_pseudo_survprob[order(pseudo$id),])		#back to original order
  return(out)
}
#----------------------------------------------------------------------------
#- A function predicts the marinal survival probability
#----------------------------------------------------------------------------
predictCondPseudo <- function(surv_train, cen_train, weight, pickTime,
                                   x_train, x_test){
  
  pseudoCond <- getPseudoConditional(surv_train, cen_train, weight, pickTime)
  
  tau=pseudoCond$tau
  s=c(0,tau[-length(tau)])
  x_train.all <- do.call(rbind, replicate(length(s), x_train, simplify=FALSE)) # repeat x_train length(s) times
  
  x_train.all <- cbind.data.frame(x_train.all, sl=rep(s,each=length(surv_train)))
  
  y_train.all <- c(pseudoCond$conditional_pseudost)
  
  dd=cbind.data.frame(y_train.all,t=rep(surv_train,length(s)), x_train.all)
  dd=dd[dd$t>dd$sl,]
  x_train.all=as.matrix(dd[,-c(1,2)])
  y_train.all=dd[,1]
  
  model <- pseudoDNN.train(x_train.all, y_train.all)
  
  
  x_test.all <- lapply(1:length(s),function(i) cbind(x_test, rep(s[i], nrow(x_test))))
  x_test.all <- Reduce(rbind, x_test.all)
  
  ypred.all <- pseudoDNN.predict(model, x_test.all)
  ypred.all <- matrix(ypred.all, nrow=nrow(x_test))
  ypred.all <- lapply(1:length(s), function(i) apply(ypred.all[,1:i, drop=FALSE], 1, prod))
  ypred.all <- Reduce(cbind, ypred.all)
  ypred.all
}


################################################
# simulations
###############################################

ptm <- proc.time()
BB=100

n=4000
p=1


pickTime=c(0.7,1.7,3.2,5.3,8.3)
res.kw=res.k=NULL

for(kk in 1:BB){
  
  set.seed(342*kk+kk+4565)
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
  
  # non-IPCW calculation
  ypred <- predictCondPseudo_ipcw(surv_train, cen_train,weight=1, pickTime, xs_train, xs_test) 
  
  res.k[[kk]]=evaluateCondPseudo(1-ypred, surv_train, cen_train, pickTime,
                                 surv_test, cen_test)
  
  # IPCW calculation
  weight=compute_weight(surv_train,cen_train,x_train)
  ypredw <- predictCondPseudo(surv_train, cen_train, weight=weight, pickTime, xs_train, xs_test) 
  
  res.kw[[kk]]=evaluateCondPseudo(1-ypredw, surv_train, cen_train, pickTime,
                                  surv_test, cen_test)
 
  print(kk)
  
}
proc.time()-ptm

