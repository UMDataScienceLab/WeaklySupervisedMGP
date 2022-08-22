mixtureCorrGPPriors <-
  function(N, nout, labs, obsLabs, ...){
    argument = list(...)
    if (!is.null(argument$obsprob)) {
      # need to add for the case that prior probs are provided
    } else {
      priors = matrix(rep(rep(1/nout, nout), N), nrow=N)
      for (n in 1:N){
        if (bdata$observedLab[n]){
          probs = rep(1e-5,nout)
          probs[labs[n]] = 1 - 1e-5*(nout-1)
          priors[n,] = probs
        } 
      }
    }
    
    return(priors)
  }