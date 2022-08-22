jitChol <-
  function ( M, Num=10, silent=FALSE ) {
    jitter <- 0
    jitter1 <- abs(mean(diag(M)))*1e-6
    eyeM <- diag( 1, nrow=length(M[,1]), ncol=length(M[1,]) )
    
    for ( i in 1:Num ) {
      ## clear the last error message
      try(stop(""),TRUE)
      
      Ch <- try( chol( M + jitter*eyeM ), silent=TRUE )
      
      nPos <- grep("not positive definite",  geterrmessage())
      
      if ( length(nPos) != 0 ) {
        jitter1 <- jitter1*10
        jitter <- jitter1
        
        if (! silent) {
          warnmsg <- paste("Matrix is not positive definite, adding",
                           signif(jitter,digits=4), "jitter!")
          warning(warnmsg)
        }
      }
      else break
    }
    
    return (list(chol=Ch, jitter=jitter))
  }
