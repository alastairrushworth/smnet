
get_optimal_smooth<-function(P.list, X.spam, XTX.spam, X.list, 
                             response, control,
                             net = F, n.linear = n.linear,  
                             sm.names, lin.names, lin.means, 
                             method, verbose, 
                             Pwee.list = Pwee.list, fixed.df){
 
  
  
  # SET UP COMPONENTS AND SIMPLE SUMMARIES REQUIRED FOR GETTING OPTIMAL PARAMETERS
  # ------------------------------------------------------------------------------
  # number of smoothing components (this is usually the same as or less than the number of smoothing pars)
  n.smooth.comps<-length(P.list)
  # vector where each element describes the number of smooth pars associated with each smooth component
  stretch_index<-vector("numeric", length = n.smooth.comps)
  for(i in 1:length(P.list)) stretch_index[i]<-ifelse(class(P.list[[i]]) == "list", length(P.list[[i]]), 1)
  # unnest the list of penalty matrices
  # P.flat    <- make_flat(P.list)
  P.flat    <- unlist(P.list, use.names = F)
  # Pwee.flat <- make_flat(Pwee.list)
  Pwee.list   <- unlist(Pwee.list, use.names = F)
  # the total number of smoothing parameters (=sum(stretch_index))
  n.smooth.pars      <-length(P.flat)
  # length of the response vector
  n         <-length(response)
  # number of columns associated with each component of the data matrix
  X.dim     <- lapply(X.list, ncol)
  # number of components in the model - includes intercept and linear variables
  n.terms   <- length(X.list)
  # total columns of the data matrix
  np        <- sum(unlist(X.dim))
  # column numbers and cumulative column numbers in the data matrix
  inds      <- unlist(X.dim)
  cum.inds  <- cumsum(inds) 
  
  # IDENTIFIABILITY CONSTRAINTS USED IN ESTIMATION OF UNIVARIATE SMOOTHS
  # ------------------------------------------------------------------------------
  blockzero  <- lapply(X.dim, make_sparse)
  idList     <- lapply(X.list, crossprodspam)
  # zero out the identifiability constraint associated with any linear components, 
  #  and include those associated with standard spline smooths
  ridgeM     <- list()
  for(z in 1:(n.linear + 1)) idList[[z]] <- blockzero[[z]] 
  for(term in 1:(n.terms - net)){
    ridgeM[[term]]<- get_block_penalty(P = idList[[term]], blockzero = blockzero, i = term)*(10^-6)
  }
  ridgeM       <- Reduce("+", ridgeM)
  # P.flat       <- c(P.flat, ridgeM)


  # CHOLESKY FACTORISATION, Xy.
  # ------------------------------------------------------------------------------
  P          <- P.flat
  Xy         <- t(X.spam) %*% response
  cholFactor <- chol.spam(XTX.spam + Reduce("+", P.flat) + ridgeM, pivot = "MMD", eps = 10^-8)
    
  # WORK OUT THE INITIAL VALUES OF THE SMOOTHING PARAMETERS FOR THE OPTIMISER
  # ------------------------------------------------------------------------------
  if(is.numeric(control$start)){
    start.vals <- control$start
  } else if(is.null(control$start)){
    start.vals <- rep(1, n.smooth.pars)
  }
  
  
  if(method %in% c("AICC", "GCV", "AIC")){
      if(is.numeric(control$approx)){
        set.seed(1)
        w <- matrix(sample(c(-1,1), n*control$approx, replace = T), nrow=n, ncol=control$approx)
        Xw<-t(X.spam) %*% w
        objective<-function(rhoArgs){
          get_crit_approx(rhoArgs, X = X.spam, 
                          XTX = XTX.spam, P = P, response = response,
                          Xy = Xy, Xw = Xw, cholFactor = cholFactor, n=n,
                          n.sm = n.smooth.pars, crit = method, identifyBit = ridgeM)
        }
      } else {
        Xw<-t(X.spam)
        objective<-function(rhoArgs){
          get_crit_exact(rhoArgs, X = X.spam, 
                         XTX = XTX.spam, P = P, response = response,
                         Xy = Xy, Xw = Xw, cholFactor = cholFactor, n=n,
                         n.sm = n.smooth.pars, crit = method, identifyBit = ridgeM)
        }
      }
      if(length(P) == 1){
        a1       <- proc.time()
        optimal  <- optimize(f = objective, upper = 30, lower = -30, tol = control$tol)
        a2       <- proc.time()
        rhoOptim <- optimal$minimum
      }
      if(length(P) > 1){
        a1       <- proc.time()
        if( ((length(P) > 2) & (!is.null(fixed.df))) | ((length(P) > 1) & (is.null(fixed.df)))){
          # get the optimal smoothing parameters
          optimal  <- optim(par = start.vals, fn = objective, method = control$optim.method, 
                            control = list(trace = control$trace, maxit = control$maxit,
                                           reltol = control$tol))
          rhoOptim <- optimal$par
        } else {
          rhoOptim <- start.vals
        }
        # if fixed.df is specified (only happens if a network smooth is specified)
        # then search for smoothing parameter values that satisfy this df level 
        # given previous optimal values for other smooth terms
          if(!is.null(fixed.df)){
            approxTF <- is.numeric(control$approx)
            optimal  <- get_fixed_df(P.flat, rhoInput = rhoOptim, 
                                     ridgeM, XTX.spam, cholFactor, info, 
                                     Xy, Xw = Xw, y = response, X.spam, 
                                     approx = approxTF,
                                     X.list, fixed.df, maxit = control$maxit)
            rhoOptim[(length(rhoOptim) - 1):length(rhoOptim)] <- optimal$par
          }
        a2       <- proc.time()
      }
    ntime   <- round(abs(a1-a2)[3], 3)
    nits    <- ifelse(length(P) > 1, optimal$counts[1], 1)
  }

  # FINALLY FIT THE MODEL USING THE ESTIMATED SMOOTHING PARAMETERS
  # ------------------------------------------------------------------------------
  for(j in 1:n.smooth.pars) P[[j]]<-P[[j]]*exp(rhoOptim[j])
  Psum     <- Reduce("+", P)
  info     <- XTX.spam + Psum + ridgeM
  U        <- update.spam.chol.NgPeyton(cholFactor, info)
  beta_hat <- backsolve.spam(U, forwardsolve.spam(U, Xy)) 
  # adjust intercept for mean centering of other covariates, if they exist.
  if(n.linear > 0) beta_hat[1] <- beta_hat[1] - sum(lin.means*beta_hat[2:(1 + n.linear)])
  fit      <- X.spam %*% beta_hat
  # get total effective numbers of degrees of freedom
  # depends on whether the approximation is in use or not
  pdof     <- vector("numeric", length = n.smooth.comps)
  if(is.numeric(control$approx)){
    left1 <- forwardsolve.spam(U, Xw) 
    ED1   <- ifelse(!is.null(dim(left1)), sum(rowMeans(left1*left1)), sum(left1*left1))
    left2 <- backsolve.spam(U, left1)
    left3 <- X.spam %*% left2
    ED2   <- sum(rowMeans(left3*left3))
    for(i in 1:n.smooth.comps){
      j        <- which(!inds == 1)[i]
      Xw_zero  <- Xw
      zero_out <- which(!(1:nrow(Xw)) %in% (cum.inds[j-1]+1):(cum.inds[j]))
      Xw_zero[zero_out,] <- 0
      right1    <- forwardsolve.spam(U, Xw_zero) 
      pdof[i]   <-  ifelse(!is.null(dim(left1)), sum(rowMeans(right1 * left1)), sum((right1 * left1)))
    }  
  } else {
    Xw    <- t(X.spam)
    left1 <- forwardsolve.spam(U, Xw) 
    ED1   <- sum(left1*left1)
    left2 <- backsolve.spam(U, left1)
    left3 <- X.spam %*% left2
    ED2   <- sum(left3*left3)
    for(i in 1:n.smooth.comps){
      j        <- which(!inds == 1)[i]
      Xw_zero  <- Xw
      zero_out <- which(!(1:nrow(Xw)) %in% (cum.inds[j-1]+1):(cum.inds[j]))
      Xw_zero[zero_out,] <- 0
      right1    <- forwardsolve.spam(U, Xw_zero) 
      pdof[i]   <- sum(right1 * left1)
    }
  }
  
  # VARIANCE AND DIAGNOSTICS FOR OUTPUT
  # ------------------------------------------------------------------------------
  sigma.sq    <- sum((response - fit)^2)/(n - 2*ED1 + ED2) # model variance
  dof.smooth  <- rep(pdof, stretch_index)
  diagnostics <- list(ED = 2*ED1 - ED2, sigma.sq = sigma.sq, fit = fit)
  
  # IF VERBOSE IS SET TO TRUE, PRINT THE CONVERGENCE STATUS FROM OPTIMISER
  # ------------------------------------------------------------------------------
  if(length(P.flat) > 1){
    if(control$verbose){
      if(optimal$convergence == 0){
        nits <- ifelse(length(P) == 1, optimal$counts, optimal$counts[1])
        cat(paste("Convergence reached in", nits, "iterations after", ntime, "s \n\n"))
      } else { cat("Warning, convergence was not reached \n\n")}   
    }
  } else if(length(P.flat) == 1) cat(paste("Convergence reached after", ntime, "s \n\n"))

  
  # PRINT THE SMOOTHING PARAMETERS, AND ASSOCIATED EFFECTIVE DIMENSIONS
  # ------------------------------------------------------------------------------
  if(!is.null(fixed.df)){
    if(method == "AICC") out <- log(sum((response - fit)^2)) + (2*(ED1+1)/(n-ED1-2))
    if(method == "AIC")  out <- log(sum((response - fit)^2)) + 2*ED1/n
    if(method == "GCV")  out <- sum((response - fit)^2)/(1-(ED1/n))^2
  }
  
  get.inner.length  <- function(L) ifelse(!class(L) == "list", 1, length(L))
  inner.length      <- lapply(P.list, get.inner.length)
  tablecol          <- max(unlist(inner.length))
  tablerow          <- length(sm.names) + net
  retPrint          <- ret <- matrix(as.factor("---"), nrow = tablerow, ncol = tablecol)
  increment         <- 1
  for(i in 1:length(P.list)){
    for(j in 1:inner.length[[i]]){
      retPrint[i,j]  <- ret[i,j] <- round(rhoOptim[increment], 2)
      increment      <- increment + 1
    }
  }
  smpar.names<-vector(length = tablecol)
  for(i in 1:tablecol) smpar.names[i]<-paste("par_", i, sep = "")
  retPrint            <- cbind(retPrint, round(unlist(pdof), 2))
  ret                 <- cbind(ret, round(unlist(pdof), 2))
  rownames(retPrint)  <- c(if(net){c(sm.names, "Network")}else{sm.names})
  rownames(ret)       <- c(if(net){c(sm.names, "Network")}else{sm.names})
  colnames(retPrint)  <- c(smpar.names, "df")
  colnames(ret)       <- c(smpar.names, "df")
  if(control$verbose){
    cat("Estimated smoothing parameters (log scale) \n")
    cat("------------------------------------------ \n")
    print(as.data.frame(retPrint)) 
  }
  function_min <- ifelse(length(P.flat) == 1, round(optimal$objective, 3), round(optimal$value, 3))
  cat(paste("\n\n", method, " = ", function_min, sep = ""))
  # RETURN MODEL OUTPUT
  # ------------------------------------------------------------------------------
  list(pars = exp(rhoOptim), beta_hat = beta_hat, ED = 2*ED1 - ED2, fit = fit, 
       ret = ret, retPrint = retPrint, sigma.sq = diagnostics$sigma.sq, U = U, lin.means = lin.means)
}


