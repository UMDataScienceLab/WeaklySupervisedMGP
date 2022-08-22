ggKernCompute <- 
  function(H, x, x2){
    
    if ( nargs()<3 ) {
      x2 = x
    }
    
    Pr = H[1]             #precision U
    Pqr = H[2]            #precision G
    sensitivity = H[3]
    Prinv = 1/Pr
    Pqinv = 1/Pqr
    Pinv = Prinv + 2*Pqinv
    
    P = 1/Pinv
    
    n2 = dist2(x, x2)
    Kbase = exp(-0.5*P*n2)
    
    K = (sensitivity^2) * Kbase
    
    return(list(K=K, Kbase=Kbase, Prinv=Prinv, Pqinv=Pqinv, P=P, n2=n2))
    
  }