mixtureCorrGPComputeVarDistMat <-
  function(model){
    M = model$nout
    
    model$DBlock = lapply(1:M, function(m){checkKernelSymmetry(model$DBlock[[m]])})
    DBlockinvRes = lapply(1:M, function(m){jitCholInv(model$DBlock[[m]])})
    model$DBlockinv = lapply(1:M, function(m){DBlockinvRes[[m]]$invM})
    
    model$varDist$varSiginv = lapply(1:M, function(m){model$DBlockinv[[m]] + diag(model$varDist$varPi[,m]*model$beta[m])})
    model$varDist$varSiginv = lapply(1:M, function(m){checkKernelSymmetry(model$varDist$varSiginv[[m]])})
    model$varDist$varSig = lapply(1:M, function(m){jitCholInv(model$varDist$varSiginv[[m]])$invM})  
    
    model$DBlockinvKfuKuuinv = lapply(1:M, function(m){model$DBlockinv[[m]]%*%t(model$KuuinvKuf[[m]])})
    
    model$varDist$varA = lapply(1:M, function(m){model$varDist$varSig[[m]] %*% model$DBlockinvKfuKuuinv[[m]]})
    model$varDist$varb = lapply(1:M, function(m){model$varDist$varSig[[m]] %*% model$By[[m]]})
    
    model$sumJSigvarSigJ = Reduce('+', lapply(1:M, function(m){model$KuuinvKuf[[m]]%*%model$DBlockinvKfuKuuinv[[m]] - t(model$DBlockinvKfuKuuinv[[m]])%*%model$varDist$varA[[m]]}))
    model$varDist$varKuuinv = model$sumJSigvarSigJ + model$Kuuinv
    model$varDist$varKuu = jitCholInv(model$varDist$varKuuinv)$invM
    
    model$sumByA = Reduce('+', lapply(1:M, function(m){t(model$By[[m]]) %*% model$varDist$varA[[m]]}) )  
    model$varDist$varmuu = model$varDist$varKuu%*%t(model$sumByA) 
    
    model$AvarKuu = lapply(1:M, function(m){model$varDist$varA[[m]] %*% model$varDist$varKuu})
    
    return(model)
  }