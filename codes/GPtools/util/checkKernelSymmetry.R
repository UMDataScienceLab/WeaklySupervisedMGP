checkKernelSymmetry <-
  function(K){
    asim.K = max(max(K-t(K)))
    if (asim.K != 0){
      K = (K + t(K))/2
    }
    
    return(K)
  }
