#----------------------------------------------------------------------------
# It trains the data in training dataset and 
# predicts marinal survival probabilities at pickTime for subjects in test dataset
# surv_train is survival time in training data
# cen_train is censoring indicator in training data
# weight is the weight in computing the pseudo probabilities 
# x_train is the covariate values in training dataset
# x_test is the covaraite values in test dataset
#----------------------------------------------------------------------------
predictSurvprob <- function(surv_train, cen_train, weight, pickTime, x_train, x_test){  
  
  pseudoCond  = getPseudoConditional(surv_train, cen_train, weight, pickTime)
  
  x_train.all <- x_train[pseudoCond$id,]
  x_train.all <- cbind(x_train.all, pseudoCond$s)
  
  y_train.all <- pseudoCond$pseudost
  
  model <- pseudoDNN.train(x_train.all, y_train.all)
  
  s <- c(0, pickTime)
  s <- s[-length(s)]
  x_test.all <- lapply(1:length(s),function(i) cbind(x_test, rep(s[i], nrow(x_test))))
  x_test.all <- Reduce(rbind, x_test.all)
  
  ypred.all <- pseudoDNN.predict(model, x_test.all)
  ypred.all <- matrix(ypred.all, nrow=nrow(x_test))
  ypred.all <- lapply(1:length(s), function(i) apply(ypred.all[,1:i, drop=FALSE], 1, prod))
  ypred.all <- Reduce(cbind, ypred.all)
  ypred.all
}
