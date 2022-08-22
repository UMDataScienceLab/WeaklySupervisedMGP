diagKernGradient <-
  function(model, partialMat, factor){
    
    # gradient with respect to variance
    g = sum(factor * diag(partialMat))
    
    return(g)
  }