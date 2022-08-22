complexLog <-
  function (x) {
    if ( is.double(x) & x>0 ) {
      y <- log(x)
    } else {
      if ( is.double(x) & x<0 )
        warning("Log of negative real number, using complex log!")
      y <- log(x+0i)
    }
    return ( y )
  }