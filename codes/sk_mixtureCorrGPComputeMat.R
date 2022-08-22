mixtureCorrGPComputeMat <-
  function(model){
    # paramsTr = spMultigpParamTrans(model)
    # betaTr = expTransform(model$beta,'atox')
    
    # compute Kuu
    Kuu.param = model$param[grep('latent', model$paramName)]
    u = model$inducingPt
    model$Kuu = rbfKernCompute(c(Kuu.param[1],1), u)$K + Kuu.param[2]*diag(model$numInducing)
    
    # compute Kfu & Kff (blockdiag)
    precisionG = model$param[grep('precisionG', model$paramName)]
    sensitivity = model$param[grep('sensitivity', model$paramName)]
    precisionU = model$param[grep('precisionU', model$paramName)]
    model$Kfu = lapply(1:model$nout, function(m){ggXGaussianKernelCompute(H=c(precisionG[m],precisionU,sensitivity[m]),
                                                                          x=model$X[[m]],
                                                                          x2=u)$K})
    model$Kff = lapply(1:model$nout, function(m){ggXggKernCompute(H=c(precisionG[m],precisionG[m],precisionU,sensitivity[m],sensitivity[m]),
                                                                  x=model$X[[m]],
                                                                  x2=model$X[[m]])$K + 1e-9*diag(length(model$X[[m]]))}) ########
    
    
    # compute Kuuinv & log(det(Kuu))
    KuuinvRes = jitCholInv(model$Kuu)
    model$Kuuinv = KuuinvRes$invM; model$sqrtKuu = KuuinvRes$chol
    model$logDetKuu = 2*sum(log(diag(model$sqrtKuu)))
    model$logDetKuuT = model$logDetKuu
    
    # compute KuuinvKuf & KfuKuuinvKuf
    model$KuuinvKuf = lapply(1:model$nout, function(m){model$Kuuinv %*% t(model$Kfu[[m]])})
    model$KfuKuuinvKuf = lapply(1:model$nout, function(m){model$Kfu[[m]] %*% model$KuuinvKuf[[m]]})
    
    # compute D (blockdiag) $ log(det(D))
    # model$Binv = lapply(1:model$nout, function(m){(1/model$beta[m])*diag(1/model$varDist$varPi[,m])})
    model$By = lapply(1:model$nout, function(m){model$beta[m] * (model$varDist$varPi[,m]) * model$m[[m]] })
    model$DBlock = lapply(1:model$nout, function(m){model$Kff[[m]] - model$KfuKuuinvKuf[[m]]})
    model$D = lapply(1:model$nout, function(m){model$Kff[[m]] - model$KfuKuuinvKuf[[m]] + (1/model$beta[m])*diag(1/model$varDist$varPi[,m])})
    model$D = lapply(1:model$nout, function(m){checkKernelSymmetry(model$D[[m]])})
    DinvRes = lapply(1:model$nout, function(m){jitCholInv(model$D[[m]])})
    model$Dinv = lapply(1:model$nout, function(m){DinvRes[[m]]$invM})
    model$sqrtD = lapply(1:model$nout, function(m){DinvRes[[m]]$chol})
    model$logDetD = lapply(1:model$nout, function(m){2*sum(log(diag(model$sqrtD[[m]])))})
    model$logDetD = Reduce('+', model$logDetD)
    model$logDetDT = model$logDetD
    
    # compute KufDinv
    model$KufDinv = lapply(1:model$nout,function(m){t(model$Kfu[[m]])%*%model$Dinv[[m]]})
    
    # compute traceDinvyy
    model$traceDinvyy = lapply(1:model$nout, function(m){sum(model$Dinv[[m]]%*%model$m[[m]]*model$m[[m]])})
    model$traceDinvyy = Reduce('+', model$traceDinvyy)
    
    # compute KufDinvy
    model$KufDinvy = Reduce('+',lapply(1:model$nout, function(m){model$KufDinv[[m]]%*%model$m[[m]]}))
    
    # compute KufDinvKfu
    model$KufDinvKfu = Reduce('+',lapply(1:model$nout, function(m){model$KufDinv[[m]]%*%model$Kfu[[m]]}))
    
    # compute A 
    model$A = model$Kuu + model$KufDinvKfu
    AinvRes = jitCholInv(model$A)
    model$Ainv = AinvRes$invM; model$sqrtA = AinvRes$chol
    model$sqrtAinv = jitChol(model$Ainv)$chol
    model$logDetA = 2*sum(log(diag(model$sqrtA))) 
    model$sqrtAinvKufDinvy = model$sqrtAinv%*%model$KufDinvy
    model$AinvKufDinvy = model$Ainv%*%model$KufDinvy
    
    return(model)
  }