mixtureCorrGPObjective <-
  function(vars, model, learn='learnhyp'){
    
    paramsTr = vars[1:model$nparam] 
    if (model$constBeta){
      betaTr = vars[(model$nparam+1):(model$nparam+1)]
    } else {
      betaTr = vars[(model$nparam+1):(model$nparam+model$nout)]
    }
    
    param = spMultigpParamTrans(model, param=paramsTr, trans='atox')
    beta = expTransform(betaTr,'atox')
    if (model$constBeta){
      beta = rep(beta, model$nout)
    }
    
    if (learn=='learnall'){
      if (model$constBeta){
        varLogPiTr = vars[(model$nparam+2):(model$nparam+1+(model$nout-1)*model$N)]
      } else {
        varLogPiTr = vars[(model$nparam+model$nout+1):(model$nparam+model$nout+(model$nout-1)*model$N)]
      }
      varLogPiTr = matrix(c(rep(0, model$N), varLogPiTr), ncol=model$nout)
      model$varDist$varPi = exp(varLogPiTr)/rowSums(exp(varLogPiTr))
      model$varDist$varLogPi = log(model$varDist$varPi)
    }
    
    model = mixtureCorrGPUpdateParam(model=model, param=param, beta=beta, updateVar=FALSE)
    
    return(-mixtureCorrGPCorLB(model=model, learn=learn)$ll)
  }