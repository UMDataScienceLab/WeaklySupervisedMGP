mixtureCorrGPCorLB <- 
  function(model, learn='learnhyp'){
    N = Reduce('+',lapply(1:model$nout, function(m){length(model$X[[m]])}))
    ll = N*log(2*pi)
    ll = ll - model$logDetKuu
    ll = ll + model$logDetA
    ll = ll + model$logDetD
    ll = ll + model$traceDinvyy
    ll = ll - t(model$KufDinvy)%*%model$Ainv%*%model$KufDinvy
    ll = -0.5*ll
    t1 = as.numeric(ll)
    # logPi = do.call('cbind', lapply(1:model$nout, function(m){rep(log(model$pi[m]), model$N)}))
    
    
    logPi = log(model$pi)
    varLogPi = model$varDist$varLogPi
    varLogPi = varLogPi - rowMaxs(varLogPi, value=T); varLogPi = varLogPi - log(rowSums(exp(varLogPi))); varPi = exp(varLogPi)
    
    
    idx.obsT = which(model$obsLab); idx.obsF = which(!model$obsLab)
    KLZ = sum(varPi[idx.obsT,]*(varLogPi[idx.obsT,]-logPi[idx.obsT,]))
    
    # need revise for more than 3 group cases
    abar = varPi[idx.obsF,]+model$Dir.a
    logB.a0 = length(idx.obsF)*log(beta(model$Dir.a, model$Dir.a))
    logB.abar = sum(log(beta(abar[,1],abar[,2])))
    KLDir = sum(varPi[idx.obsF,]*(varLogPi[idx.obsF,])) - logB.abar + logB.a0
    
    if (learn=='learnall' && model$constPriorProb) KLZ = 0 
    t3 = Reduce("+",lapply(1:model$nout, function(m){sum( 0.5*((1-varPi[,m])*log(2*pi*1/model$beta[m]) - varLogPi[,m]) )} ) )
    # t3 = 0 ########
    ll = ll - KLZ - KLDir + t3
    # ll = t3 ####
    # ll = varLogPi[1,1] ###########
    return(list(ll=as.numeric(ll), t1=t1, KLZ=KLZ, KLDir=KLDir, t3=t3))
  }