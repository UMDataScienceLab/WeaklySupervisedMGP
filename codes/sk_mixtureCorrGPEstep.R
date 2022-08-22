# Update elements 1:n in qZ at iteration n, 
# then 100 (can be altered by maxiter) iterations more updating all the elements in qZ

mixtureCorrGPEstep <-
  function(model, maxiter=100){
    
    maxiteration = maxiter + model$N   
    
    varLogPi = model$varDist$varLogPi
    varLogPi = varLogPi - rowMaxs(varLogPi, value=T); varLogPi = varLogPi - log(rowSums(exp(varLogPi))); varPi = exp(varLogPi)
    logPi = do.call('cbind', lapply(1:model$nout, function(m){rep(log(model$pi[m]), model$N)}))
    
    old = -Inf; rNA=rep(NA, maxiteration); corLBTab=data.frame(ll=rNA, KLZ=rNA, t3=rNA)
    for (iter in 1:maxiteration){
      
      # update qz
      varmu = lapply(1:model$nout, function(m){model$varDist$varb[[m]] + model$varDist$varA[[m]]%*%model$varDist$varmuu})
      varVar = lapply(1:model$nout, function(m){unlist(lapply(1:model$N,function(n){model$AvarKuu[[m]][n,] %*%t(model$varDist$varA[[m]])[,n] })) + diag(model$varDist$varSig[[m]])})
      PI = lapply(1:model$nout, function(m){(model$m[[m]] - varmu[[m]])^2 + varVar[[m]]}) 
      tmp = lapply(1:model$nout, function(m){ -0.5*model$beta[m]*PI[[m]] -0.5*log(2*pi*(1/model$beta[m])) + logPi[,m] })
      tmp = do.call('cbind', tmp)
      varLogPi[1:min(iter, model$N),] = tmp[1:min(iter, model$N),]
      varLogPi = varLogPi - rowMaxs(varLogPi, value=T); varLogPi = varLogPi - log(rowSums(exp(varLogPi))); varPi = exp(varLogPi)
      
      for (m in 1:model$nout){
        varPi[,m][varPi[,m]>0.99] = rep(0.99, length(varPi[,m][varPi[,m]>0.99] )) 
        varPi[,m][varPi[,m]<0.01] = rep(0.01, length(varPi[,m][varPi[,m]<0.01] )) 
      }
      varLogPi = log(varPi)
      
      # rowMaxs(varPi, value=T) > 1-(1e-)
      model$varDist$varLogPi = varLogPi
      model$varDist$varPi = varPi
      
      # update qfqu
      model = mixtureCorrGPComputeMat(model)
      model = mixtureCorrGPComputeVarDistMat(model)
      
      corLBResult = mixtureCorrGPCorLB(model)
      corLBTab$ll[iter] = corLBResult$ll; corLBTab$KLZ[iter] = corLBResult$KLZ; corLBTab$t3[iter] = corLBResult$t3; 
      corLB = corLBResult$ll; model$correctedLB = corLB
      
      print(paste0('iteration: ', iter, ' corrected LB: ', corLB))
      
      if ( iter > model$N & abs(corLB - old) < abs(corLB)*1e-8 ){
        break
      } else {
        old = corLB
      }
    }
    
    return(model)
  
  }