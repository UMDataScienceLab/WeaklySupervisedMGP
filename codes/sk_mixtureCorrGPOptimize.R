mixtureCorrGPOptimize <-
  function(model, EM=TRUE, iter=300, optimizer='SCG'){
    
    # Unexponentiating the parameters
    paramsTr = spMultigpParamTrans(model, param=model$param, trans='xtoa')
    betaTr = expTransform(model$beta,'xtoa')
    if (model$constBeta){
      betaTr = betaTr[1]
    } 
    
    # vars: variables to be optimized
    vars = c(paramsTr, betaTr) # before transforming (for unconstrained feasible region)
    
    # Check gradient before optimize
    # gradCheck(vars=vars, func=mixtureCorrGPObjective, grad=mixtureCorrGPGradient, model=model)
    
    optimOptions = list(optimizer=optimizer,
                        maxit=iter,
                        display=1,
                        gradCheck=TRUE,
                        xtol=1e-4,
                        fnTol=1e-4)
    
    if (EM){
      ##### EM algorithm: optimizing both variational and original distribution
      print("****** Start EM algorithm ******")
      old = -Inf
      for (i in 1:6) {
        
        
        ## E-step: update varational distributions
        model = mixtureCorrGPEstep(model, maxiter=100)
        
        ## M-step: update hyperparameters
        # Unexponentiating the parameters
        paramsTr = spMultigpParamTrans(model, param=model$param, trans='xtoa')
        betaTr = expTransform(model$beta,'xtoa')
        if (model$constBeta){
          betaTr = betaTr[1]
        } 
        # vars: variables to be optimized
        vars = c(paramsTr, betaTr) # before transforming (for unconstrained feasible region)
        model = mixtureCorrGPMstep(model, vars, options=optimOptions)
        
        if (abs(model$correctedLB - old) < abs(model$correctedLB) * 1e-4) break
        else old = model$correctedLB
      }
    }
    
    ### final optimization in terms of every variables
    paramsTr = spMultigpParamTrans(model, param=model$param, trans='xtoa')
    betaTr = expTransform(model$beta,'xtoa')
    if (model$constBeta){
      betaTr = betaTr[1]
    } 
    varqLogPi = as.vector((model$varDist$varLogPi - model$varDist$varLogPi[,1])[,(2:model$nout)])
    # vars: variables to be optimized
    vars = c(paramsTr, betaTr, varqLogPi) # before transforming (for unconstrained feasible region)
    
    if (model$constPriorProb){
      varplogPi = (log(model$pi) - log(model$pi)[1])[2:model$nout]
      vars = c(vars, varplogPi)
    }
    
    
    # gradCheck(vars=vars, func=mixtureCorrGPObjective, grad=mixtureCorrGPGradient, model=model, learn='learnall')
    # stop()
    model = mixtureCorrGPMstep(model, vars, options=optimOptions, learn='learnall')
    
    return(model)
  }