#' @importFrom spam cbind.spam

predictSSNobject <- function(object, newdata){
  
  #   put all of the components of new design matrix in a list
  newX <- vector("list")
  
  # extract the prediction locations and the prediction data matrix
  if(is.null(newdata)){
    # rids     <- object$ssn.object@data$rid[object$ssn.object@data$netID == object$internals$netID]
    # rid_ord  <- rids[order(rids)]
    dfPred   <- getSSNdata.frame(object$ssn.object, "preds")
    dfPred   <- dfPred[dfPred$netID == object$internals$netID, ]
    ridPred  <- as.numeric(as.character(dfPred[,"rid"]))
  } else {
    dfPred   <- newdata
    dfPred   <- dfPred[dfPred$netID == object$internals$netID, ]
    ridPred  <- as.numeric(as.character(dfPred[,"rid"]))
  }
  
  # get the linear components if it exists
  if(object$internals$n.linear > 0){
    linVarbs  <- dfPred[,colnames(object$internals$variables)[1:object$internals$n.linear]]
    newX      <- c(newX,  list(linear = cbind(1, linVarbs)))
  } else {
    newX      <- c(newX, list(linear = rep(1, length(ridPred))))
  }
  
  #   get the smooth additive components if there are any  
  if(object$internals$n.smooth > 0){
    smTerms      <- object$internals$sm.terms.names
    smDesign     <- vector("list")
    dataOriginal <- data.frame(object$internals$variables)
    for(i in 1:length(smTerms)){
      oldVariable  <- dataOriginal[smTerms[[i]]]
      new_bit      <- construct_bspline_basis(as.matrix(dfPred[smTerms[[i]]]), 
                                                     dimensions = object$internals$sm.basis[[i]], 
                                                     cyclic = object$internals$sm.cyclic[[i]], 
                                                     range.variables = oldVariable)
      newX         <- c(newX, new_bit)
    }
  }
  
  
  
  if(object$internals$net){
    # construct network component
    ridPredRemap  <- re_map_rid(rid_vector = ridPred, all_rid = as.integer(object$internals$adjacency$rid_bid[, 1]))
    newX <- c(newX, 
              spam(
                x = list(i   = 1:length(ridPredRemap), 
                         j   = ridPredRemap, 
                         val = rep(1, length(ridPredRemap))), 
                nrow = length(ridPred), 
                ncol = ncol(object$internals$X.list[[length(object$internals$X.list)]])))
  }
  
  # put together the linear, smooth and network components
  newX           <- lapply(newX, as.matrix)
  Xstar          <- Reduce("cbind.spam", newX)
  predictions    <- as.numeric(Xstar %*% object$internals$beta_hat)
  left1          <- forwardsolve.spam(object$internals$U, t(object$internals$X.spam))
  left2          <- backsolve.spam(object$internals$U, left1)
  vec            <- Xstar %*% left2
  predictions.se <- sqrt((1 + rowSums(vec*vec))*object$internals$sigma.sq)
  list(predictions = predictions, predictions.se = predictions.se)   
  }
