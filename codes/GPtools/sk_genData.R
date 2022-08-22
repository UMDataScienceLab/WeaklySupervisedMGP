genData <-
  function(dataset, ..., seedVal){
    
    if ( nargs()<2 ) {
      seedVal = 1e+6;
    } 
    set.seed(seedVal)
    print(sprintf("Seed is set to %f", seedVal))
    
    # setting for parameters to generate data
    if (dataset=='ggToy'){
      
      args <- list(...)
      nout = args$nout
      nlf = args$nlf
      precisionG = args$precisionG
      sensitivity = args$sensitivity
      precisionU = args$precisionU
      N = args$N
      if ( is.null(args$bias) ) {
        bias = rep(0, nout)
      } else {
        bias = args$bias  
      }
      
      if (nout != length(precisionG)) {stop('Wrong input')}
      if (length(precisionG) != length(sensitivity)) {stop('Wrong input')}
      if (length(precisionU) > 1) {stop('Not implemented')}
      
      x = seq(-1,1,length.out=N)
      X = list()
      for (i in 1:nout) X[[i]] = x
      
      mu = rep(0,N*nout)
      for (i in 1:nout){
        if (bias[i] != 0){
          mu[((i-1)*N+1):(i*N)] = rep(bias[i], N)
        }
      }
      
      fullCov = spMultiFullBlockCompute(nout=nout, X=X, precisionU=precisionU, 
                                        precisionG=precisionG, sensitivity=sensitivity)
      Ytrue = matrix(c(gaussSamp(mu=mu, Sigma=fullCov, numSamps=1)),ncol=nout)
      Ytrue = lapply(1:nout, function(i) Ytrue[,i])
      Yobs = list()
      for (i in 1:nout) Yobs[[i]] = Ytrue[[i]] + 0.1*sqrt(var(Ytrue[[i]]))*rnorm(N)
      return(list(X=X, y=Yobs, ytrue=Ytrue, bias=bias))
    }
    
    else if (dataset=='ggToyMissing'){
      bdata = genData('ggToy', 
                      nout = 2,
                      nlf = 1,
                      precisionG=c(120,200),
                      sensitivity=c(4,5),
                      precisionU=100,
                      N = 500,
                      bias = c(0,2.5),
                      seedVal=seedVal)
      arguments = list(...) # diffX = TRUE
      if (is.null(arguments$diffX)){
        diffX = FALSE
      } else diffX = arguments$diffX
      
      nout = length(bdata$X)
      gdata = list('X'=bdata$X, 'y'=bdata$y)
      
      N = 120
      missingUnit = c(1)
      missingLb = 0; missingUb = 0
      # missingIdx = c(1:length(bdata$X[[1]]))[bdata$X[[i]] < missingLb | bdata$X[[i]] > missingUb]
      
      obsIdx = lapply(1:nout, function(n){sort(sample(c(1:length(bdata$X[[n]])), size=N))})
      obsBool = lapply(1:nout, function(n){ c(1:length(bdata$X[[n]])) %in% obsIdx[[n]] })
      obsMissBool = obsBool
      
      if (!diffX) {
        obsIdx = lapply(1:nout, function(n){obsIdx[[1]]})
        obsBool = lapply(1:nout, function(n){obsBool[[1]]})
        obsMissBool = obsBool
      }
      
      # obsIdx = sort(sample(c(1:length(bdata$X[[1]])), size=N))
      # obsBool = c(1:length(bdata$X[[1]])) %in% obsIdx
      # obsMissBool = obsBool
      for (n in 1:nout){
        for (i in 1:length(bdata$X[[n]])){
          if ( missingLb < bdata$X[[n]][i] & bdata$X[[n]][i] < missingUb) obsMissBool[[n]][i] = FALSE
        }
      }
      
      for (i in 1:nout){
        if (i %in% missingUnit) {
          gdata$X[[i]] = gdata$X[[i]][ obsMissBool[[i]] ]
          gdata$y[[i]] = gdata$y[[i]][ obsMissBool[[i]] ]
        }
        else {
          gdata$X[[i]] = gdata$X[[i]][obsBool[[i]]]
          gdata$y[[i]] = gdata$y[[i]][obsBool[[i]]]
        }
      }
      
      gdata$trueVal=bdata
      
      return(gdata)
    }
    
    else if (dataset=='ggToyMissingNoLab'){
      bdata = genData('ggToyMissing',..., seedVal=seedVal)
      arguments = list(...)
      sparseOut = c(1) #
      sparseFactor = c(1)
      if (!is.null(arguments$sparseOut)) sparseOut = arguments$sparseOut
      if (!is.null(arguments$sparseFactor)) sparseFactor = arguments$sparseFactor
      nout = length(bdata$X)
      
      sparseList=list()
      for (i in 1:length(sparseOut)){
        nX = length(bdata$X[[sparseOut[i]]])
        sparseListTmp = rep(FALSE, nX)
        numTrue = as.integer(nX*sparseFactor[[i]])
        sparseListTmp[1:numTrue] = rep(TRUE, numTrue)
        permSparseListTmp = sample(sparseListTmp)
        bdata$X[[sparseOut[i]]] = bdata$X[[sparseOut[i]]][permSparseListTmp]
        bdata$y[[sparseOut[i]]] = bdata$y[[sparseOut[i]]][permSparseListTmp]
      }
      gdata = list()
      gdata$X = as.vector(unlist(bdata$X))
      gdata$y = as.vector(unlist(bdata$y)) 
      labList = list()
      for (i in 1:nout){
        labList[[i]] = rep(i, length(bdata$X[[i]]))
      }
      gdata$lab = unlist(labList)
      gdata$nout = nout
      gdata$trueVal = bdata$trueVal
      
      return(gdata)
    }
    
    else if (dataset=='ggToyMissingPartialLab'){
      gdata = genData('ggToyMissingNoLab',..., seedVal=seedVal)
      nout = length(gdata$trueVal$X)
      arguments = list(...)
      obsRangeLB = min(gdata$X); obsRangeUB = max(gdata$X) # need a revision for the different ranges over groups   
      obsGroup = c(1:nout)
      obsRate = rep(0,nout)
      if (!is.null(arguments$obsRange)) {
        obsRangeLB=rep(0, nout); obsRangeUB=rep(0, nout)
        for (i in 1:nout){
          obsRangeLB[i] = arguments$obsRange[[i]][1]; obsRangeUB[i] = arguments$obsRange[[i]][2]
        }
      }
      if (!is.null(arguments$obsGroup)) obsGroup = arguments$obsGroup
      if (!is.null(arguments$obsRate)) obsRate = arguments$obsRate
      
      gdata$observedLab = rep(FALSE,length(gdata$X))
      for (i in 1:nout){
        X.i.idx = which(gdata$lab==i)
        X.i = gdata$X[X.i.idx]
        X.i.idx.inRange = X.i.idx[ X.i>=obsRangeLB[i] & X.i<=obsRangeUB[i] ] 
        nX.i.inRange = length(X.i.idx.inRange)
        nObs = floor(nX.i.inRange*obsRate[i])
        X.i.obsidx = sample(X.i.idx.inRange, size=nObs)
        gdata$observedLab[X.i.obsidx] = TRUE
      }
      return(gdata)
    }
    
  }