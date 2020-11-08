#' Specify Network Smoother in Formulae
#' 
#' @description This function specifies all of the information required to smooth parameters 
#' over the segments of a stream network using an adjacency matrix, and a vector of flow weights.  
#' 
#' @param adjacency A sparse adjacency matrix of class \code{spam} that describes the flow 
#' connectedness of the stream network.  \code{adjacency} is typically obtained from a call 
#' to \link{get_adjacency}
#' @param weight A character string indicating the column name of a numeric vector of flow 
#' weights contained in the \code{data.object} that has been passed to \link{smnet}.  
#' Defaults to "\code{autoShreve}" which automatically constructs a weighting based on 
#' Shreve order, useful if data does not include an appropriate weight.  For more information
#'  on choosing appropriate \code{weight} inputs from a given data set, see \link{show_weights}.
#' @param fixed.df Positive scalar indicating a fixed number of degrees of freedom to allocate 
#' to the stream network component, overriding the criterion minimisation for this component. 
#' Under the default setting, NULL, the degrees of freedom are chosen automatically.
#' 
#' @return A list combining the processed input components above.  
#' For internal use within \code{smnet}.
#' \itemize{
#' \item{\code{adjacency}: Sparse adjacency matrix}
#' \item{\code{weight}: Numeric vector of flow weights}
#' \item{\code{netID}: Integer identifying network of interest}
#' }
#' 
#' @author Alastair Rushworth
#' @export
network <- function(adjacency = NULL,  weight = "autoShreve", fixed.df = NULL){
  ret <- list(adjacency = adjacency,  weight = weight, fixed.df = fixed.df, by = "NA")
  class(ret) <- "network.spec"
  ret
}
