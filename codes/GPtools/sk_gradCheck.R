gradCheck <-
  function(vars, func, grad, model, ...){
    
    eps = 1e-6
    
    tmpModel = model
    # func = spMultigpObjective ### should be removed
    # grad = spMultigpGradient ###
     
    nparams = length(vars)
    deltaf = rep(0, nparams)
    step = rep(0, nparams)
    
    for (i in 1:(nparams)){
      step[i] = 1.0
      fplus = func(vars=vars+(eps*step), model=tmpModel, ...) # fplus = func(vars=vars+(eps*step), model=tmpModel, learn='learnall')
      fminus = func(vars=vars-(eps*step), model=tmpModel, ...) # fminus = func(vars=vars-(eps*step), model=tmpModel, learn='learnall')
      deltaf[i] = 0.5*(fplus-fminus)/eps
      step[i] = 0.0
    }
    gradient = grad(vars=vars, model=tmpModel, ...)
    delta = deltaf - gradient
    
    if (model$constBeta){
      paramName = c(model$paramName, sprintf("Beta %d",c(1:1)))
    } else paramName = c(model$paramName, sprintf("Beta %d",c(1:model$nout)))
    
    paramName = c(paramName, sprintf("Additional var %d", c(1:(nparams - length(paramName))) ))
    
    for (i in 1:nparams){
      print(sprintf("param: %-25s analytic: %-12f numerical: %-12f diff: %-12f", paramName[i], gradient[i], deltaf[i], delta[i]))
    }
  }