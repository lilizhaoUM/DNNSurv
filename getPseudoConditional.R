#-----------------------------------------------------------------------------------------------------------------------------------
# get pseudo conditional survival probabilities 
# t is the survival time
# d is the censoring indicator
# qt is a vector of time points that are used to divide the time interval
# output has subject id (id) and time points (s) and pseudo conditional survival probabilities (pseudost) for subject=id and at time s
#------------------------------------------------------------------------------------------------------------------------------------
getPseudoConditional <- function(t, d, qt){
  #browser()
  s <- c(0, qt)  
  n=length(t)
  ns=length(s)-1  # the number of intervals
  D <- do.call(cbind, lapply(1:ns, function(j)  (s[j] < t) * (t <= s[j+1]) * (d == 1)))
  R <- do.call(cbind, lapply(1:ns, function(j) ifelse(s[j] < t, 1, 0)))
  Delta<-do.call(cbind, lapply(1:ns, function(j) pmin(t,s[j+1])-s[j]))
  
  # format into long formate
  dd.tmp=cbind.data.frame(id=rep(1:n,ns),s=rep(c(0,qt[-length(qt)]), each=n), y=c(R*Delta),d=c(D))
  
  dd=dd.tmp[dd.tmp$y>0,]
  pseudost=rep(NA, nrow(dd))
  for (j in 1:ns){
    index= (dd$s==s[j])
    dds=dd[index,]
    pseudost[index]=pseudosurv(time=dds$y, event=dds$d, tmax=s[j+1]-s[j])$pseudo
    print(j)
  }
  dd$pseudost=pseudost  
  
  return(dd[,c(1,2,5)])
}



