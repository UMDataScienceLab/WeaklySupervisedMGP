mixtureCorrGPMstep <-
  function(model, vars, options, learn='learnhyp'){
    
    optimizer = options$optimizer
    
    # if (learn=='learnhyp'){
    if (optimizer=='SCG'){
      optimRes = SCG(x=vars, fn=mixtureCorrGPObjective, grad=mixtureCorrGPGradient, options=options, model=model, learn=learn)
      optimVars = optimRes$xmin
    }
    else if (optimizer=='nloptr'){
      library('nloptr')
      nloptrOpt <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = options$maxit, print_level=2)
      optimRes = nloptr(x0=vars, eval_f=mixtureCorrGPObjective, eval_grad_f=mixtureCorrGPGradient, opts=nloptrOpt, model=model, learn=learn)
      optimVars = optimRes$solution
    }
    # }
    
    optimParamTr = optimVars[1:model$nparam]
    if (model$constBeta) {
      optimBetaTr = optimVars[(model$nparam+1):(model$nparam+1)]
      optimBetaTr = rep(optimBetaTr, model$nout)
    } else optimBetaTr = optimVars[(model$nparam+1):(model$nparam+model$nout)]
    
    optimParam = spMultigpParamTrans(model=model, param=optimParamTr, trans='atox')
    optimBeta = expTransform(optimBetaTr,'atox')
    
    if (learn=='learnall'){
      if (model$constBeta){
        varLogPiTr = optimVars[(model$nparam+2):(model$nparam+1+(model$nout-1)*model$N)]
      } else {
        varLogPiTr = optimVars[(model$nparam+model$nout+1):(model$nparam+model$nout+(model$nout-1)*model$N)]
      }
      varLogPiTr = matrix(c(rep(0, model$N), varLogPiTr), ncol=model$nout)
      model$varDist$varPi = exp(varLogPiTr)/rowSums(exp(varLogPiTr))
      model$varDist$varLogPi = log(model$varDist$varPi)
    }
    
    model = mixtureCorrGPUpdateParam(model, optimParam, optimBeta)
    model$correctedLB = -optimRes$objective
    
    return(model)
    
  }