mixtureCorrGPCreate <-
  function(X, y, options){
  model = list()
  
  ## input data
  model$nout = options$nout
  model$nlf = options$nlf       # number of latent function
  model$X = list()   # X data
  model$y = list()   # y data
  for (i in 1:model$nout){
    model$X[[i]] = X
    model$y[[i]] = y
  }
  model$N = length(X)
  minx = min(model$X[[1]]); maxx = max(model$X[[1]])
  model$obsLab = options$obsLab
  
  ## parameters
  model$fixInducing = options$fixInducing 
  model$constPriorProb = options$constPriorProb
  model$numInducing = options$numActive
  
  if (model$fixInducing){
    model$nparam = 2*model$nout*model$nlf + 2*model$nlf
    model$param = rep(1, model$nparam)
  } else {
    model$nparam = 2*model$nout*model$nlf + 2*model$nlf + model$numInducing
    model$param = rep(1, model$nparam)
    model$param[(model$nparam-model$numInducing+1):model$nparam] = seq(minx, maxx, length.out = model$numInducing)
  }
  
  model$paramName = rep(NA, model$nparam)
  model$paramName[1:(model$nout*model$nlf)] = sprintf("precisionG out %d",c(1:model$nout))
  model$paramName[(model$nout*model$nlf+1):(2*model$nout*model$nlf)] = sprintf("sensitivity out %d",c(1:model$nout))
  model$paramName[(2*model$nout*model$nlf+1):(2*model$nout*model$nlf+model$nlf)] = sprintf("precisionU latent %d",c(1:model$nlf))
  model$paramName[(2*model$nout*model$nlf+model$nlf+1):(2*model$nout*model$nlf+2*model$nlf)] = sprintf("w variance latent %d",c(1:model$nlf))
  if (!model$fixInducing){
    model$paramName[(model$nparam-model$numInducing+1):model$nparam] = sprintf("inducing %d",c(1:model$numInducing))
  } 
  factorInducing = 0.1; med = (maxx-minx)/2 
  model$inducingPt = seq(minx - factorInducing*med, maxx + factorInducing*med, length.out=model$numInducing)
  
  model$param[grep('w variance latent', model$paramName)] = exp(-2)
  
  model$constBeta = options$constBeta
  if (model$constBeta && length(options$beta)==1){
    model$beta = rep(options$beta, model$nout)
  } else {
    model$beta = options$beta
  }
  
  model$betaTrans = options$betaTrans
  
  model$nonnegative = rep(FALSE, model$nparam)
  model$nonnegative[grep("precision", model$paramName)] = TRUE
  model$nonnegative[grep("variance", model$paramName)] = TRUE
  
  if (is.null(options$bias)) {
    model$bias = rep(0,model$nout)
  } else model$bias = options$bias
  
  if (is.null(options$mu)) {
    model$mu = lapply(1:model$nout, function(m){rep(0,length(model$X[[m]]))})
  } else model$mu = options$mu
  
  # if (is.null(options$pi)) {
  #   model$pi = rep(1/model$nout, model$nout)
  # } else model$pi = options$priorProb
  
  if (is.null(options$priorProb)) {
    model$pi = matrix(rep(rep(1/model$nout, model$nout), model$N), nrow=model$N)
  } else model$pi = options$priorProb
  
  model$varDist = list()
  
  if (!is.null(options$varPi)){
    model$varDist$varPi = options$varPi
  } else {
    varPi = matrix(runif(model$nout*model$N)+10, ncol=model$nout)
    model$varDist$varPi = varPi/rowSums(varPi)
  }
  model$varDist$varLogPi = log(model$varDist$varPi)
  
  model$Dir.a = options$Dir.a
  
  ## calculate important matrix
  model = spMultigpComputeM(model)
  model = mixtureCorrGPComputeMat(model)
  model = mixtureCorrGPComputeVarDistMat(model)
  
  model$VarDistMatUpdated = TRUE
  
  model$correctedLB = NA
  
  return(model)
}