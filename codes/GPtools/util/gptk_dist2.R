dist2 <-
  function (x, x2) {
    xdim <- dim(as.matrix(x))
    x2dim <- dim(as.matrix(x2))
    
    xMat <- array(apply(as.matrix(x*x),1,sum), c(xdim[1], x2dim[1]))
    x2Mat <- t(array(apply(as.matrix(x2*x2),1,sum), c(x2dim[1], xdim[1])))
    
    if ( xdim[2] != x2dim[2] )
      stop("Data dimensions are not matched.")
    
    n2 <-   xMat+x2Mat-2*tcrossprod(x, x2)
    
    return (n2)
  }