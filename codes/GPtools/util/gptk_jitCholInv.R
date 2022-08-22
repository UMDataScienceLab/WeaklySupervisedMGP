jitCholInv <-
  function ( M, Num=15, silent=FALSE ) {
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
    
    invCh <- try (solve( Ch, eyeM ), silent=TRUE)
    
    if ( class(invCh) == "try-error" ) {
      return (NaN)
    }
    else {
      invM <- invCh %*% t(invCh)
      
      if ( jitter == 0 ) {
        ans <- list(invM=invM, jitter=jitter, chol=Ch)
      }
      else ans <- list(invM=invM, jitM=M+jitter*eyeM , jitter=jitter, chol=Ch)
      
      return (ans)
    }
  }