ggXggKernGradient <-
  function(H, x, x2, partialMat){
    
    # compute kernel
    kernRes = ggXggKernCompute(H=H, x=x, x2=x2)
    K = kernRes$K
    Kbase = kernRes$Kbase
    Pqrinv = kernRes$Pqrinv
    Psrinv = kernRes$Psrinv
    Prinv = kernRes$Prinv
    P=kernRes$P
    fsensq = kernRes$fsensq
    fsenss = kernRes$fsenss
    dist = kernRes$n2
    
    preFactorPqr = 1/(2*Pqrinv + Prinv)
    preFactorPsr = 1/(2*Psrinv + Prinv)
    preFactorPr = (1/(2*Pqrinv + Prinv)) + (1/(2*Psrinv + Prinv))
    
    dim = matrix(1, nrow=length(x), ncol=length(x2)) 
    
    matGradPr = sum( 0.5*K*(Prinv*Prinv*(dim*P - 0.5*dim*preFactorPr - P*dist*P))*partialMat )
    matGradPqr = sum( 0.5*K*(Pqrinv*Pqrinv*(dim*P - dim*preFactorPqr - P*dist*P))*partialMat )
    matGradPsr = sum( 0.5*K*(Psrinv*Psrinv*(dim*P - dim*preFactorPsr - P*dist*P))*partialMat )
    GradSensq = fsensq * sum( Kbase*partialMat )
    GradSenss = fsenss * sum( Kbase*partialMat )
    
    return(list(matGradPr=matGradPr,
                matGradPqr=matGradPqr,
                matGradPsr=matGradPsr,
                GradSensq=GradSensq,
                GradSenss=GradSenss))
    
  }