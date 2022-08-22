ggXGaussianKernGradient <-
  function(H, x, x2, partialMat){
    
    kernRes = ggXGaussianKernelCompute(H=H, x=x, x2=x2)
    K=kernRes$K
    Kbase=kernRes$Kbase
    Pqrinv=kernRes$Pqrinv
    Prinv=kernRes$Prinv
    P=kernRes$P
    fsens=kernRes$fsens
    dist=kernRes$n2
    
    dim = matrix(1, nrow=length(x), ncol=length(x2))
    preFactorPqr = 1/(2*Pqrinv + Prinv)
    preFactorPr = 1/Prinv + (1/(2*Pqrinv + Prinv))
    
    matGradPr = sum( 0.5*partialMat*K*(Prinv*Prinv*(dim*P - 0.5*dim*preFactorPr- P*dist*P)) )
    matGradPqr = sum( 0.5*partialMat*K*(Pqrinv*Pqrinv*(dim*P - dim*preFactorPqr- P*dist*P)) )
    
    gradSenss = fsens*sum(partialMat*Kbase);
    
    return(list(matGradPr=matGradPr,
                matGradPqr=matGradPqr,
                gradSenss=gradSenss))
    
  }