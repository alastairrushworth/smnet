#' Summarise a Stream Network Model
#' 
#' @description Generate summaries of linear and smooth components of 
#' an \code{smnet} object.
#' 
#' @param object An object of class \code{smnet}.
#' @param verbose Logical.  Whether to print summaries to the console.
#' @param ...  other arguments passed to summary.
#' 
#' @return List object with components
#' \itemize{
#' \item{1}{\code{linear.terms}: the linear components of the fitted model}
#' \item{2}{\code{smooth.terms}: the values of the smoothing parameters on the log scale, and the partial degrees of freedom associated with each smooth component.  Note: Network components always have two smoothing parameters, where the second is a (usually small) ridge parameter.}
#' }
#' 
#' @author Alastair Rushworth
#' @seealso   \code{\link[=smnet]{smnet}}, \code{\link{plot.smnet}}
#' @export
#' @importFrom stats pt
#' @aliases summary.smnet
#' @examples   

#' # As an example, create a simulated SSN object
#' # Save the object to a temporary location
#' set.seed(12)
#' ssn_path <- paste(tempdir(), "/example_network", sep = "")
#' 
#' # If example network doesn't already exist, then attempt to create it
#' # Otherwise, read from the temporary directory
#' example_network <- try(importSSN(ssn_path, o.write = TRUE), silent = TRUE)
#' if('try-error' %in% class(example_network)){
#' example_network <- createSSN(
#'       n            = 50,
#'       obsDesign    = binomialDesign(200), 
#'       predDesign   = binomialDesign(50),
#'       importToR    = TRUE, 
#'       path         = ssn_path,
#'       treeFunction = iterativeTreeLayout
#'   )
#' }
#' 
#' ## create distance matrices, including between predicted and observed
#' createDistMat(example_network, "preds", o.write=TRUE, amongpred = TRUE)
#' 
#' ## extract the observed and predicted data frames
#' observed_data            <- getSSNdata.frame(example_network, "Obs")
#' prediction_data          <- getSSNdata.frame(example_network, "preds")
#' 
#' ## associate continuous covariates with the observation locations
#' #  data generated from a normal distribution
#' obs                      <- rnorm(200)
#' observed_data[,"X"]      <- obs
#' observed_data[,"X2"]     <- obs^2
#' 
#' ## associate continuous covariates with the prediction locations
#' #  data generated from a normal distribution
#' pred                     <- rnorm(50) 
#' prediction_data[,"X"]    <- pred
#' prediction_data[,"X2"]   <- pred^2
#' 
#' ## simulate some Gaussian data that follows a 'tail-up' spatial process
#' sims <- SimulateOnSSN(
#'   ssn.object      = example_network, 
#'   ObsSimDF        = observed_data, 
#'   PredSimDF       = prediction_data,  
#'   PredID          = "preds",  
#'   formula         = ~ 1 + X,
#'   coefficients    = c(1, 10),
#'   CorModels       = c("Exponential.tailup"), 
#'   use.nugget      = TRUE,
#'   CorParms        = c(10, 5, 0.1),
#'   addfunccol      = "addfunccol")$ssn.object
#' 
#' 
#' ## extract the observed and predicted data frames, now with simulated values
#' sim1DFpred         <- getSSNdata.frame(sims, "preds")
#' sim1preds          <- sim1DFpred[,"Sim_Values"]
#' sim1DFpred[,"Sim_Values"] <- NA
#' sims               <- putSSNdata.frame(sim1DFpred, sims, "preds")
#' 
#' # create the adjacency matrix for use with smnet
#' adjacency    <- get_adjacency(
#'   ssn_path, 
#'   net = 1
#' )
#' 
#' # not run - plot the adjacency matrix
#' # display(adjacency[[1]])
#' 
#' # sometimes it is useful to see which variables are valid network weights 
#' # in the data contained within the SSN object
#' show_weights(sims, adjacency)
#' 
#' # fit a penalised spatial model to the stream network data
#' # Sim_Values are quadratic in the X covariate.  To highlight 
#' # the fitting of smooth terms, this is treated as non-linear 
#' # and unknown using m().
#' mod_smn       <- smnet(formula = Sim_Values ~ m(X) + m(X2) + 
#'                        network(adjacency = adjacency, weight = "shreve"), 
#'                        data.object = sims, netID = 1)
#' 
#' # not run - plot different summaries of the model
#' plot(mod_smn, type = "network-covariates")
#' plot(mod_smn, type = "network-segments", weight = 4, shadow = 2)
#' plot(mod_smn, type = "network-full", weight = 4, shadow = 2)
#' 
#' # obtain predictions at the prediction locations and plot 
#' # against true values
#' preds <- predict(mod_smn, newdata = getSSNdata.frame(sims, "preds"))
#' plot(preds$predictions, sim1preds)
#' 
#' # obtain summary of the fitted model
#' summary(mod_smn)
#' 
#' #' # delete the simulated data
#' unlink(ssn_path, recursive = TRUE)


summary.smnet<-function(object, verbose = TRUE, ...){
  adjacency  <- object$internals$adjacency
  n.smooth   <- object$internals$n.smooth
  retPrint   <- object$internals$retPrint
  n.linear   <- object$internals$n.linear
  X.list     <- object$internals$X.list
  beta_hat   <- object$internals$beta_hat
  U          <- object$internals$U
  sigma.sq   <- object$internals$sigma.sq
  XTX.spam   <- object$internals$XTX.spam
  X.spam     <- object$internals$X.spam
  lin.names  <- object$internals$lin.names
  ED         <- object$internals$ED
  fit        <- object$internals$fit
  n          <- length(object$internals$response)
  response   <- object$internals$response

  # Summarise the linear part of the model
  linearExists <- n.linear + 1
  if(linearExists > 0){
    getcol    <- function(M) ifelse(ncol(M) == "NULL", 1, ncol(M)) 
    cov.dims  <- lapply(X.list, getcol)    
    inds      <- unlist(cov.dims)
    cum.inds  <- cumsum(inds)
    n.cov     <- length(X.list)
    unit.locations <- which(inds == 1)
    
    # get standard errors 
    left1         <- forwardsolve.spam(U, t(X.spam))
    left2         <- backsolve.spam(U, left1)
    separs        <- sqrt(sigma.sq*rowSums(left2^2)[unit.locations])
    pars          <- beta_hat[unit.locations]
    ret.mat       <- cbind(pars, separs, pars/separs, pt(abs(pars/separs), df = n - ED, lower.tail = FALSE) * 2)
    ret           <- as.data.frame(matrix(prettyNum(ret.mat, digits = 5), nrow = length(pars), byrow = F))
    
    rownames(ret) <- c("(Intercept)", lin.names) 
    colnames(ret) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    if(verbose){
      cat("\n-----------------------------------------------\nLinear terms:\n-----------------------------------------------\n")
      print(ret)
    }
  }
  if(verbose){
    if(n.smooth + !is.null(adjacency) > 0){
      cat("\n\n-----------------------------------------------\nSmooth terms:\n-----------------------------------------------\n")
      print(as.data.frame(retPrint))
      cat("\n")
    }
  }
  
  # return the summarise quietly
  invisible(list(linear.terms = ret, smooth.terms = as.data.frame(retPrint)))
}




  