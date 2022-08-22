mixtureCorrGPUpdateParam <-
  function(model=model, param=param, beta=beta, updateVar=TRUE){
    
    model$param = param
    model$beta = beta
    model = mixtureCorrGPComputeMat(model)
    model$VarDistMatUpdated = FALSE
    
    if (updateVar == TRUE){
      model = mixtureCorrGPComputeVarDistMat(model)
      model$VarDistMatUpdated = TRUE
    }
    
    return(model)
  }