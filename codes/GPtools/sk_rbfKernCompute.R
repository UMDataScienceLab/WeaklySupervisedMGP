rbfKernCompute <-
  function (H, x, x2=NULL) {
    
    precisionU = H[1]
    variance = H[2]
    
    if ( nargs() < 3 ) {
      n2 <- dist2(x,x)
    } else {
      n2 <- dist2(x,x2)
    }
    
    wi2 <- 0.5*precisionU
    Kbase <- exp(-n2*wi2)
    k <- variance*Kbase
    
    return (list(K=k, Kbase=Kbase, n2=n2))
  }
