#' @importFrom spam update.spam.chol.NgPeyton
#' @importFrom spam backsolve.spam
#' @importFrom spam forwardsolve.spam
#' @importFrom stats optim

get_fixed_df <- function(P.flat, rhoInput, approx, ridgeM, XTX.spam, cholFactor, info, Xy, Xw, X.spam, X.list, y, maxit, fixed.df){
  X.dim     <- lapply(X.list, ncol)
  n.terms   <- length(X.list)
  np        <- sum(unlist(X.dim))
  inds      <- unlist(X.dim)
  cum.inds  <- cumsum(inds) 
  P         <- P.flat
  
  get_df <- function(rho, rhoInput, ridgeM, Xw, X.spam, cholFactor, approx, XTX.spam, P, Xy, n.terms, cum.inds){
    rhoInput[(length(rhoInput) - 1):length(rhoInput)] <- rho
    for(j in 1:length(rhoInput)) P[[j]]<-P[[j]]*exp(rhoInput[j])
    info     <- XTX.spam + Reduce("+", P) + ridgeM
    U        <- tryCatch(update.spam.chol.NgPeyton(cholFactor, info), error=function(e) e, warning=function(w) w)
    if("warning" %in% class(U) | any(abs(rhoInput) > 20)){
      out <- rnorm(1, mean = 10^20)
    } else {
      beta_hat <- backsolve.spam(U, forwardsolve.spam(U, Xy))  
      left1    <- forwardsolve.spam(U, Xw) 
      Xw_zero  <- Xw
      zero_out <- which(!(1:np) %in% (cum.inds[n.terms-1]+1):(cum.inds[n.terms]))
      Xw_zero[zero_out,] <- 0
      right1    <- forwardsolve.spam(U, Xw_zero) 
      if(approx){
        out <- abs(ifelse(!is.null(dim(left1)), sum(rowMeans(right1 * left1)), sum((right1 * left1))) - fixed.df)
      } else {
        out <- abs(sum(right1 * left1) - fixed.df)
      }
      
    }
    out
  }
  
  opt_df <- optim(par = c(5, -5), fn = get_df, rhoInput = rhoInput, P = P, ridgeM = ridgeM, XTX.spam = XTX.spam, 
                  cholFactor = cholFactor, approx = approx, Xy = Xy, Xw = Xw, n.terms = n.terms, cum.inds = cum.inds, X.spam = X.spam, 
                  control = list(maxit = maxit))
  # need to work out which pars 
  opt_df
}
