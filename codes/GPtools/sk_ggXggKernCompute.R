ggXggKernCompute <-
  function(H, x, x2){
    
    if ( nargs()<3 ) {
      x2 = x
    }
    
    Pqr = H[1]  # Precision G
    Psr = H[2]  # Precision G
    Pr = H[3]   # Precision U
    Pqrinv = 1/Pqr
    Psrinv = 1/Psr
    Prinv = 1/Pr
    Pinv = Pqrinv + Psrinv + Prinv
    P = 1/Pinv
    
    Sqr = H[4]  # Sensitivity G
    Ssr = H[5]  # Sensitivity G
    
    n2 = dist2(x, x2)
    fNumPqr = (2*Pqrinv + Prinv)^(1/4)
    fNumPsr = (2*Psrinv + Prinv)^(1/4)
    fDen = prod(Pinv)^(1/2)
    factor = Sqr * Ssr * fNumPqr * fNumPsr / fDen
    fsensq = Ssr * fNumPqr * fNumPsr / fDen
    fsenss = Sqr * fNumPqr * fNumPsr / fDen
    
    Kbase = exp(-0.5 * P * n2)
    K = factor * Kbase
    
    return(list(K=K, Kbase=Kbase, Pqrinv=Pqrinv, Psrinv=Psrinv, Prinv=Prinv, P=P,
                fsensq=fsensq, fsenss=fsenss, n2=n2))
    
  }