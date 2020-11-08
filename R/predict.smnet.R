#' Make Predictions Using a Fitted Stream Network Model.
#' 
#' @description Get predictions and standard errors for a new set of spatial locations 
#' and associated covariate values from a model fitted by \code{smnet()}.
#' 
#' @param object Object of class \code{smnet}, usually the result of a call to \code{smnet()}
#' @param newdata \code{data.frame} containing new covariate values at which to make predictions.
#' @param ... other arguments passed to \code{predict.smnet()}
#' 
#' @return List object with components
#' \itemize{
#' \item{\code{predictions}: vector of predictions corresponding to prediction points in 
#' the original \code{SpatialStreamNetwork} input object}
#' \item{\code{predictions.se}: vector of prediction standard errors}
#' }
#' 
#' @author Alastair Rushworth
#' @export
#' 
predict.smnet <- function(object, newdata = NULL, ...){
  out <- predictSSNobject(object = object, newdata = newdata)
  return(out)
}