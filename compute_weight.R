#---------------------------------------------------------------------------------------------------
# -compute weights using Cox regression for the ipcw pseudo appraoch 
# -output is n*ns matrix: n is the number of subjects and ns is the number of unique event times
#--------------------------------------------------------------------------------------------------
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
