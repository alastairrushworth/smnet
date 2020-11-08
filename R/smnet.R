#' Additive Modelling for Stream Networks
#' 
#' @description Fits (Gaussian) additive models to river network data based on the 
#' flexible modelling framework described in O'Donnell et al. (2014).  Data must be 
#' in the form of a \code{SpatialStreamNetwork} object as used by the \code{SSN} 
#' package (Ver Hoef et al., 2012).   Smoothness of covariate effects is represented
#'  via a penalised B-spline basis (P-splines) and parameter estimates are obtained 
#'  using penalised least-squares.  Optimal smoothness is achieved by optimization 
#'  of \code{AIC}, \code{GCV} or \code{AICC}.  The formula interpreter used for 
#'  penalised additive components is modelled on the code found in the package 
#'  \code{mgcv}. 
#' 
#' @details 
#' \code{control} is a list whose elements control smoothness selection: 
#' \code{maxit} limits the number of iterations made by the optimiser (default  = 500).  
#' \code{approx}, positive integer, if specified indicates the number of Monte-Carlo 
#' samples to collect using an approximate version of performance criterion when direct 
#' evaluation is slow - this may be much faster if the network has a large number of 
#' segments or the data is large, for example \code{approx  = 100} often works well 
#' (defaults to \code{NULL}).  \code{checks}, logical, specifies whether additivity 
#' checks should be performed on the input weights default = \code{TRUE}.  
#' \code{trace}, default = 0, if set to 1, \code{optim} will print progress.   
#' \code{tol}, the relative tolerance of \code{optim}.  \code{optim.method}, 
#' the optimiser - default is "Nelder-Mead" see \code{?optim} for details.
#' 
#' @param formula A formula statement similar to those used by \code{lm} and \code{mgcv:::gam}.  
#' Smooth functions are represented with \code{m(..., k = 20)} function, up to 2d smooth 
#' interactions are currently supported.  Spatial network components are specified using 
#' \code{network(...)} function.  Further details can be found in \link{m} and \link{network} 
#' and the examples below.
#' @param data.object An  object of class \code{SpatialStreamNetwork}.
#' @param netID Integer indicating the particular stream network to model, generally only 
#' user-specified when multiple networks are contained within \code{data.object}, default is 1. 
#' @param method Character string determining the performance criterion for choosing optimal 
#' smoothness, options are \code{"AICC"}, \code{"AIC"} or \code{"GCV"}. 
#' @param control A list of options that control smoothness selection via optimisation.  
#' See 'Details'.
#' 
#' @return Object of class \code{smnet} with components
#' \itemize{
#' \item{fitted.values}{vector of fitted values}
#' \item{residuals}{vector of residuals: response minus fitted values}
#' \item{coefficients}{vector of parameter estimates}
#' \item{R2}{R^2 statistic}
#' \item{R2.adj}{adjusted R^2 statistic}
#' \item{df.residual}{residual degrees of freedom}
#' \item{ssn.object}{unchanged SSN input data object}
#' \item{internals}{model objects for internal use by other functions}
#' }
#' 
#' @author Alastair Rushworth
#' @seealso   \code{\link[=smnet]{get_adjacency}}, \code{\link{plot.smnet}}
#' @references 
#' Ver Hoef, J.M.,  Peterson, E.E., Clifford, D., Shah, R.  (2012)   SSN: An R Package for Spatial Statistical Modeling on Stream Networks 
#' O' Donnell, D., Rushworth, A.M., Bowman, A.W., Scott, E.M., Hallard, M.  (2014) Flexible regression models over river networks.  Journal of the Royal Statistical Society: Series C (Applied Statistics). 63(1) 47--63.
#' Reinhard Furrer, Stephan R. Sain (2010). spam: A Sparse Matrix R
#' Package with Emphasis on MCMC Methods for Gaussian Markov Random
#' Fields. Journal of Statistical Software, 36(10), 1-25. URL: http://www.jstatsoft.org/v36/i10/
#' @import SSN
#' @export
#' @aliases smnet
#' @keywords network
#' @keywords p-spline
#' @examples   
#' # As an example, create a simulated SSN object
#' # Save the object to a temporary location
#' set.seed(12)
#' ssn_path <- paste(tempdir(), "/example_network", sep = "")
#' 
#' # If example network doesn't already exist, then attempt to create it
#' # Otherwise, read from the temporary directory
#' example_network <- try(importSSN(ssn_path, predpts = 'preds', o.write = TRUE), silent = TRUE)
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
#
#' 
#' # plot the simulated network structure with prediction locations
#' # plot(example_network, bty = "n", xlab = "x-coord", ylab? = "y-coord")
#' 
#' ## create distance matrices, including between predicted and observed
#' createDistMat(example_network, "preds", o.write = TRUE, amongpred = TRUE)
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
#' # delete the simulated data
#' unlink(ssn_path, recursive = TRUE)




smnet <- function(formula, data.object, netID = 1, method = "AICC", control = NULL)
{  
  
  adjacency <- NULL
  default.control = list(trace = 0, maxit = 1000, checks = T, start = NULL, approx = NULL, 
                         verbose = TRUE, tol = 10^-6, optim.method = "Nelder-Mead")
  if(!is.null(control)) for(i in 1:length(control)) default.control[names(control)[i]] <- list(control[[i]])
  
  
  # ERROR CHECK AND PROCESS THE control INPUT
  if(!is.null(default.control)){
    maxit <- default.control$maxit
    if(!is.null(maxit)) if(!is.wholenumber(maxit)) stop("maxit should be a positive integer or NULL", domain = NA, call. = FALSE)
    approx <- default.control$approx
    if(!is.null(approx)) if(!is.wholenumber(approx)) stop("approx should be a positive integer or NULL", domain = NA, call. = FALSE) 
    do.checks <- default.control$checks
    if(!is.null(do.checks)) if(!is.logical(do.checks)) stop("check should be logical or NULL", domain = NA, call. = FALSE) 
    if(is.null(do.checks)) do.checks <- TRUE
  }
  
  # ERROR CHECK AND PROCESS THE method INPUT
  if(!method %in% c("AICC", "AIC", "GCV")){
    err.msg <- "method should be either AIC, AICC or GCV"
    stop(err.msg, domain = NA, call. = FALSE)
  }
  
  # ERROR CHECK AND PROCESS THE DATA INPUT
  cls.data <- class(data.object)
  if(cls.data == "SpatialStreamNetwork"){
    data            <- getSSNdata.frame(data.object, Name = "Obs")
    data_netID      <- as.numeric(as.character(data$netID))
    data            <- data[data_netID == netID, ]
  } else if(cls.data == "data.frame"){
    data <- data.object
  } else stop("input data.object must be of class SpatialStreamNetwork or data.frame", domain = NA, call. = FALSE)
  
  # ERROR CHECK AND PROCESS THE FORMULA INPUT
  # if there is an error on the formula, print and error and show the list of available data columns
  formulaout <- try(formula_stuff <- get_formula_stuff(formula = formula, data = data), silent=T)
  if(class(formulaout) == "try-error"){
    var.names <- c(t(cbind(sort(as.character(names(data))), "\n")))
    err.msg   <- paste(c(var.names, "\nFormula error: check covariates are from the above"), collapse = " ")
    stop(err.msg, domain = NA, call. = FALSE)
  } 
  
  # construct the vector of weights, apply an ordering if appropriate
  cls.sm  <- lapply(formulaout$gp$smooth.spec, class)
  net     <- "network.spec" %in% cls.sm
  fixed.df  <- NULL
  if(net){
    networkSpecNo  <- which(cls.sm == "network.spec")
    networkSpec    <- formulaout$gp$smooth.spec[[networkSpecNo]]
    adjacency      <- networkSpec$adjacency
    if(!is.null(adjacency)){
      if(class(data.object) == "SpatialStreamNetwork"){
#         netID           <- networkSpec$netID
        # get the rids out of the SSN object associated with netID
        data            <- data[as.numeric(as.character(data$netID)) == netID,]
        rid_data        <- data.object@data
        rid_data        <- rid_data[as.numeric(as.character(rid_data$netID)) == netID,]
        rid_data$rid    <- re_map_rid(rid_data$rid, as.numeric(adjacency$rid_bid[,1]))
        if(nrow(rid_data) == 0) stop("No data associated with that netID", domain = NA, call. = FALSE)
        ord             <- order(rid_data$rid)
        response.locs   <- as.numeric(as.character(data$rid))
        response.locs   <- re_map_rid(response.locs, as.numeric(adjacency$rid_bid[,1]))
        weight          <- networkSpec$weight
        # set the maximum degrees of freedom for the network smooth (default of half number sites)
        if(!is.null(networkSpec$fixed.df)) fixed.df <- networkSpec$fixed.df
        # this happens if the user requests automatic stream weighting based on Shreve order
        if(weight == "autoShreve"){
          weight       <- adjacency_to_shreve(adjacency = adjacency)
          if(default.control$checks){
            weight.type  <- check_weight(adjacency, weight)
            if(weight.type == "additive"){
              weight       <- get_shreve_weights(adjacency = adjacency$adjacency, shreve.order = as.matrix(weight))  
            } else if(weight.type == "unrecognised"){
              stop("supplied weight vector is neither additive or a network weighting", domain = NA, call. = FALSE)
            }
          }
        } else {
          # this happens if the user requests the use of a ready made stream weighting or additive function
          weight       <- try(as.vector(as.matrix(rid_data[weight])), silent = T)
          if(class(weight) == "try-error"){
            # if the requested weight does not exist, show an error
            stop("Requested weight vector not found.  \n", domain = NA, call. = FALSE)
          }
          weight       <- weight[ord]
          weight.type  <- check_weight(adjacency, weight)
          if(weight.type == "additive"){
            weight     <- get_shreve_weights(adjacency = adjacency$adjacency, shreve.order = as.matrix(weight))  
          } else if(weight.type == "unrecognised"){
            stop("supplied weight vector is neither additive or a network weighting", domain = NA, call. = FALSE)
          }
        }
      }
      if(class(data.object) == "data.frame"){
        if(is.character(networkSpec$locs)){
          response.locs <- try(data[networkSpec$locs], silent = T)
          if(response.locs == "try-error"){
            stop("No data found for locs provided")
          }
        } else if(is.numeric(networkSpec$locs)){
          response.locs <- networkSpec$locs
        } else stop("locs should be numeric vector or a character string")
        weight <- networkSpec$weight
        if(!is.numeric(weight)) stop("If no data is not an SSN object, weight must be a vector", domain = NA, call. = FALSE)
        ord <- rid_data <- netID <- NULL
      }
      # check that weights vector and adjacency are the same dimension
      if(!nrow(adjacency$adjacency) == length(weight)) stop("weights and adjacency have different dimensions")
    } else {ord <- rid_data <- netID <- NULL; stop("network smooth requested, but adjacency matrix is missing", domain = NA, call. = FALSE)}
  } else ord <- rid_data <- weight <- NULL

  # create model objects
  model_objects <- 
    get_model_objects(
      formula = formula, data = data, 
      adjacency = adjacency, 
      response.locs = response.locs, 
      weight = weight,  
      rid_data = rid_data, 
      netID = netID, ord = ord, 
      control = default.control, 
      formulaout = formulaout, 
      fixed.df = fixed.df
    )
  
  #  choose optimal smooth parameters using box constrained Nelder-Mead search
  opt <- with(model_objects, {
    if((n.smooth == 0) && (!net)){
      get_lm_fit(
        X.spam = X.spam, 
        X.list = X.list, 
        XTX.spam = XTX.spam, 
        response = response, 
        lin.means = lin.means, 
        n.linear = n.linear
      )         
    } else {
      get_optimal_smooth(
        P.list = P.list, 
        X.spam = X.spam, 
        X.list = X.list, 
        XTX.spam = XTX.spam, 
        response = response, 
        control = default.control, 
        net = net, 
        n.linear = n.linear,
        lin.names = lin.names, 
        sm.names = sm.names, 
        lin.means = lin.means,
        method = method, 
        Pwee.list = Pwee.list, 
        fixed.df = fixed.df
      )   
    }
  }
  )
  
  # calculate and print R2
  response      <- model_objects$response
  fitted.values <- opt$fit

  R2     <- 1 - sum((response - fitted.values)^2)/sum((response - mean(response))^2)
  np     <- ifelse(is.null(opt$ED), ncol(model_objects$X.spam), opt$ED)
  R2.adj <- R2 - ((1-R2)*np)/(nrow(model_objects$X.spam) - np - 1)
  if(default.control$verbose) cat(paste("   n = ", (nrow(model_objects$X.spam)), sep = ""))
  if(default.control$verbose) cat(paste("  R2.adj", " = ", round(R2.adj, 3), sep = ""))

  # create output list
  outputList               <- vector("list")
  outputList$response      <- response
  outputList$fitted.values <- fitted.values
  outputList$residuals     <- response - fitted.values
  outputList$coefficients  <- opt$beta_hat
  outputList$R2            <- R2
  outputList$R2.adj        <- R2.adj
  outputList$df.residual   <- length(response) - np
  outputList$ssn.object    <- data.object
  outputList$internals     <- c(opt, model_objects)
  class(outputList)        <- "smnet"
  invisible(outputList)
}

