#----------------------------------------------------------------------------
# A two-layer deep neural network model implemented in keras with default hyperparameters 
# hyper-parameters in this model worked very well in simulation studies in the paper
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
