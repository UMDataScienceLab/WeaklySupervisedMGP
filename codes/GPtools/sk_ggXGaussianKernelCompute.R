ggXGaussianKernelCompute <-
  function(H, x, x2){
    
    if ( nargs()<3 ) {
      x2 = x
    }
    
    Pqr = H[1]    # precision G
    Pr = H[2]     # precision U
    Pqrinv = 1/Pqr
    Prinv = 1/Pr
    Pinv = Pqrinv + Prinv
    P = 1/Pinv
    sensitivity = H[3] # sensitivity G
    
    n2 = dist2(x,x2)
    fNumPqr = (2*Pqrinv + Prinv)^(1/4)
    fNumPr = Pr^(-1/4)
    fDen = Pinv^(1/2)
    factor = sensitivity * fNumPqr * fNumPr / fDen
    fsens = fNumPqr * fNumPr / fDen
    
    Kbase = exp(-0.5*P*n2)
    K = factor * Kbase
    
    
    return(list(K=K, Kbase=Kbase, Pqrinv=Pqrinv, Prinv=Prinv, P=P, fsens=fsens, n2=n2))
  }