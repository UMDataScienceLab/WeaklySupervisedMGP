whiteKernCompute <-
  function (H, x, x2) {
    
    variance = H[1]
    
    if ( nargs()<3 ) {
      xdim <- dim(as.array(x))[1]
      k <- variance*diag(1, nrow=xdim, ncol=xdim)
    } else {
      x1dim <- dim(as.array(x))[1]
      x2dim <- dim(as.array(x2))[1]
      k <- matrix(0, nrow=x1dim, ncol=x2dim)
    }
    return (k)
  }
