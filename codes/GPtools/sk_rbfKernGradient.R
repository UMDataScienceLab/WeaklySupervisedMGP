rbfKernGradient <-
  function(H, x, x2=NULL, partialMat){
    
    kernRes = rbfKernCompute(H=H, x=x, x2=x2)
    K = kernRes$K
    n2 = kernRes$n2
    
    matGradPrecisionU = sum( -0.5*K*n2*partialMat )
    
    return(list(matGradPrecisionU=matGradPrecisionU))
  }