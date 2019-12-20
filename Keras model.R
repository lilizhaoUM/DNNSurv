#--------------------------------------------------------------------------------------
# There are two example neural network models used in the paper.
# You should tunning the hyperparamters to find the best neural network for your own study.
#-------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
# Example of a neural network with two hidden layers implemented in R keras in the paper
#-----------------------------------------------------------------------------------------
pseudoDNN.train <- function(x_train, y_train){
  
  model <- keras_model_sequential() %>%
    layer_dense(units = 8, kernel_regularizer = regularizer_l2(0.0001), activation = "tanh",
                input_shape = dim(x_train)[[2]]) %>%
    layer_dense(units = 4, kernel_regularizer = regularizer_l2(0.01),
                activation = "tanh") %>%
    layer_dense(units = 1, activation='sigmoid')
  
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

#----------------------------------------------------------------------------------------
# Another example of a neural network with one hidden layer implemented in R keras in the paper
#-----------------------------------------------------------------------------------------
pseudoDNN.train <- function(x_train, y_train){
  # use selu instead of relu for some studies
  model <- keras_model_sequential() %>%
    layer_dense(units=16,  activation = "selu",bias_initializer = initializer_constant(0.0),
                input_shape = dim(x_train)[[2]]) %>%
    layer_dropout(rate = 0.2) %>%
    layer_dense(units = 1, activation='sigmoid')
  
  model %>% compile(
    optimizer = optimizer_rmsprop(lr = 0.001),
    loss = "mse",
    metrics = c("mae")
  )
  model %>% fit(x_train, y_train,
                epochs = 1000, batch_size =256,
                verbose = 0)
  model
}

#----------------------------------------------------------------------------
#prediction based on a keras model
#----------------------------------------------------------------------------
pseudoDNN.predict <- function(model, x_test){
  ypred <- model %>% predict(x_test)
  ypred
}
