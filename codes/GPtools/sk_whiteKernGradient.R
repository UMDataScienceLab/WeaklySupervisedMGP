whiteKernGradient <-
  function(model, partialMat){
    
    # gradient with respect to variance
    
    g = sum(diag(partialMat))
    
    return(g)
  }