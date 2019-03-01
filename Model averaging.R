#===================================================================================================
# tuneCondPseudoCindexCV function tune parameters based on c index using five-fold cv.
# It selects five best sets of hyperparameters to train the model five times, output five models. 
#===================================================================================================

tuneCondPseudoCindexCV <- function(x_train, surv_train, cen_train, weight, pickTim, 
                                   tuneParas, tolerance=0.0001, maxit=20){
  
  #weight=compute_weight(surv_train,cen_train,cov.cen)
  pseudoCond <- getPseudoConditional(surv_train, cen_train, weight, pickTime)
  
  G <- 5
  n <- nrow(x_train)
  folds <- cut(seq(1,n),breaks=G,labels=FALSE)
  
  x_train.all <- y_train.all <- x_test.all <- vector("list", G)
  for(i in 1:G){
    pick.train.id <- (1:n)[folds!=i]
    pick.train.pseudoCond  <- pseudoCond[pseudoCond$id %in% pick.train.id,]
    x_train.all[[i]] <- x_train[pick.train.pseudoCond$id,]
    x_train.all[[i]] <- cbind(x_train.all[[i]], pick.train.pseudoCond$s)
    y_train.all[[i]] <- pick.train.pseudoCond$pseudost
    
    pick.test  <- (1:n)[folds==i]
    x_test <- x_train[pick.test,]
    s <- c(0, pickTime)
    s <- s[-length(s)]
    x_test.all[[i]] <- lapply(1:length(s),function(i) cbind(x_test, rep(s[i], nrow(x_test))))
    x_test.all[[i]] <- Reduce(rbind, x_test.all[[i]])
  }
  
  
  #initilize the tuning parameters
  pick.tune <- NULL
  pick.tune.new <- which(tuneParas[,1]==8 & tuneParas[,2]==4 &
                           tuneParas[,3]==0.0001 & tuneParas[,4]==0.01 &
                           tuneParas[,5]==0.0025)
  pick.tune <- c(pick.tune, pick.tune.new)
  
  pick.tune.new <- sample((1:nrow(tuneParas))[-pick.tune], 4)
  pick.tune <- c(pick.tune, pick.tune.new)
  
  pick.tune.best5 <- pick.tune
  cindex.pick.tune.best5 <- rep(0, 5)
  for(i in 1:5){
    ypred.all <- vector("list", G)
    for(j in 1:G){
      model <- pseudoDNN.train.threelayers(x_train.all[[j]], y_train.all[[j]], 
                                           tuneParas[pick.tune[i],1], tuneParas[pick.tune[i],2],
                                           tuneParas[pick.tune[i],3], tuneParas[pick.tune[i],4],
                                           tuneParas[pick.tune[i],5])
      
      
      ypred.all[[j]] <- pseudoDNN.predict(model, x_test.all[[j]])
      ypred.all[[j]] <- matrix(ypred.all[[j]], nrow=nrow(x_test.all[[j]])/length(pickTime))
      ypred.all[[j]] <- lapply(1:length(s), function(i) apply((ypred.all[[j]])[,1:i, drop=FALSE], 1, prod))
      ypred.all[[j]] <- Reduce(cbind, ypred.all[[j]])
      ypred.all[[j]] <- 1 - ypred.all[[j]]
    }
    
    ypred.all <- Reduce(rbind, ypred.all)
    cindex <- rep(0, length(pickTime))
    for(k in 1:length(pickTime)){
      cindex[k] <- concordance.index(x=ypred.all[,k], surv.time=surv_train,
                                     surv.event=cen_train, method="noether")$c.index
    }
    cindex.pick.tune.best5[i] <- mean(cindex)
  }
  
  pick.cindex <- cindex.pick.tune.best5
  pick5 <- order(cindex.pick.tune.best5, decreasing = TRUE)
  pick.tune.best5 <- pick.tune.best5[pick5]
  cindex.pick.tune.best5 <- cindex.pick.tune.best5[pick5]
  
  
  flag <- 1
  it <- 0
  while(flag){
    pick.tune.new <- sample((1:nrow(tuneParas))[-pick.tune], 1)
    pick.tune <- c(pick.tune, pick.tune.new)
    
    ypred.all <- vector("list", G)
    for(j in 1:G){
      model <- pseudoDNN.train.threelayers(x_train.all[[j]], y_train.all[[j]], 
                                           tuneParas[pick.tune.new,1], 
                                           tuneParas[pick.tune.new,2],
                                           tuneParas[pick.tune.new,3], 
                                           tuneParas[pick.tune.new,4],
                                           tuneParas[pick.tune.new,5])
      
      
      ypred.all[[j]] <- pseudoDNN.predict(model, x_test.all[[j]])
      ypred.all[[j]] <- matrix(ypred.all[[j]], nrow=nrow(x_test.all[[j]])/length(pickTime))
      ypred.all[[j]] <- lapply(1:length(s), function(i) apply((ypred.all[[j]])[,1:i, drop=FALSE], 1, prod))
      ypred.all[[j]] <- Reduce(cbind, ypred.all[[j]])
      ypred.all[[j]] <- 1 - ypred.all[[j]]
    }
    
    
    ypred.all <- Reduce(rbind, ypred.all)
    cindex <- rep(0, length(pickTime))
    for(k in 1:length(pickTime)){
      cindex[k] <- concordance.index(x=ypred.all[,k], surv.time=surv_train,
                                     surv.event=cen_train, method="noether")$c.index
    }
    print(cindex)
    
    
    if(!any(is.na(cindex))){ 
      pick.cindex <- c(pick.cindex, mean(cindex))
      
      if(abs(mean(cindex.pick.tune.best5) - mean(cindex)) < tolerance)
        flag <- 0
      
      if (mean(cindex) > cindex.pick.tune.best5[5]){
        pick.tune.best5 <- c(pick.tune.best5, pick.tune.new)
        cindex.pick.tune.best5 <- c(cindex.pick.tune.best5, mean(cindex))
        pick5 <- order(cindex.pick.tune.best5, decreasing = TRUE)[1:5]
        pick.tune.best5 <- pick.tune.best5[pick5]
        cindex.pick.tune.best5 <- cindex.pick.tune.best5[pick5]
      }
      
      it <- it+1
      if(it > maxit)
        flag <- 0
      
      print(it)
    }
    else{
      it <- it+1
      if(it > maxit)
        flag <- 0
      
      print(it)
    }
  }
  
  list(pick.tune.best5=pick.tune.best5,
       cindex.pick.tune.best5=cindex.pick.tune.best5,
       it=it, pick.tune)
  
}


#===========================================================================================
# predictCondPseudoTune function estimates the survival 
# probability by averaging the survival probabilities from the five models.
#===========================================================================================                           
predictCondPseudoTune <- function(surv_train, cen_train, weight, pickTime,
                              x_train, x_test, cov.cen, paras){
  
    pseudoCond  = getPseudoConditional(surv_train, cen_train, weight=weight, pickTime)
    
    x_train.all <- x_train[pseudoCond$id,]
    x_train.all <- cbind(x_train.all, pseudoCond$s)
    
    y_train.all <- pseudoCond$pseudost
    
    model <- pseudoDNN.train.threelayers(x_train.all, y_train.all,
                                         paras[1], paras[2],
                                         paras[3], paras[4], paras[5])
    
    s <- c(0, pickTime)
    s <- s[-length(s)]
    x_test.all <- lapply(1:length(s),function(i) cbind(x_test, rep(s[i], nrow(x_test))))
    x_test.all <- Reduce(rbind, x_test.all)
    
    ypred.all <- pseudoDNN.predict(model, x_test.all)
    ypred.all <- matrix(ypred.all, nrow=nrow(x_test))
    ypred.all <- lapply(1:length(s), function(i) apply(ypred.all[,1:i, drop=FALSE], 1, prod))
    ypred.all <- Reduce(cbind, ypred.all)
    ypred.all <- 1 - ypred.all
    
    ypred.all
}
