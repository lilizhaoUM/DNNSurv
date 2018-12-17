
#--------------------------------------------------------------------------------------------------
# get pseudo conditional survival probabilities 
# time is the survival time
# delta is the censoring indicator
# weight is the weight matrix for the IPCW pseudo probabilities. To compute non-IPCW pseduo values, set weight=1
# tau is a vector of time points that are used to divide the time interval
# output conditional_pseudost is a matrix of pseudo conditional survival probabilities with rows 
# equal to the nunmber of subjects and columns equal to the number of time points
#-------------------------------------------------------------------------------------------------
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
