get_penalty <- function(ncol, cyclic = F){
  if(!cyclic){
    zvec     <- rep(0, ncol)
    zvec[1]  <- -1; zvec[2] <- 2; zvec[3] <- -1
    circ     <- circulant.spam(zvec)
    lastrw   <- circ[nrow(circ),]
    circ     <- circ[-nrow(circ),]
    circ     <- circ[-nrow(circ),]
    out      <- t(circ) %*% circ
  } else {
    # varb_support <- abs(diff(range(varb)))
    # knotspace    <- varb_support/(ncol + 6 + 1)                  
    out          <- crossprod(circulant.spam(c(2* cos((2 * pi) / ncol), -1, rep(0, (ncol - 3)), -1)))
  }
  out
}
