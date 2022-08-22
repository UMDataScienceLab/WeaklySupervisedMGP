ggKernGradient <-
  function(H, x, x2=NULL, partialMat){
    
    if ( nargs()<4 ) {
      x2 = x
    }
    
    kernRes = ggKernCompute(H, x=x, x2=x2)
    K = kernRes$K
    Kbase=kernRes$Kbase
    Prinv=kernRes$Prinv
    Pqinv=kernRes$Pqinv
    P=kernRes$P
    dist=kernRes$n2
    sensitivity = H[3]
    
    matGradPqr = - sum( partialMat * K * (Pqinv*P)^2 * dist )
    matGradPr  = - 0.5 * sum( partialMat * K * ((Prinv*P)^2*dist) );
    gradSensq = 2 * sensitivity * sum(Kbase * partialMat)
    
    return(list(matGradPqr=matGradPqr,
                matGradPr=matGradPr,
                gradSensq=gradSensq))
    
  }