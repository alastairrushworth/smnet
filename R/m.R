#' Specify Smooth Terms in Formulae
#' 
#' @description Function used to set up univariate or bivariate smooth terms 
#' based on P-splines, for use within a call to \code{\link{smnet}}.
#' 
#' @param ... one or more variables for creating P-spline smooths.
#' @param k integer defining the number of uniformly spaced B-spline basis functions 
#' for the smooth, default is 10.  For 2d (and higher) smooths, this is the marginal 
#' basis size.
#' @param cyclic logical vector indicating whether the smooth should be cyclic.  
#' Based on the harmonic smoother of Eiler and Marx (2004)
#' 
#' @return List object with components
#' \itemize{
#' \item{\code{term}: character vector of the names of the variables involved in the smooth to be set up}
#' \item{\code{bs.dim}: number of B-spline basis functions to be used in the smooth}
#' }
#' 
#' @references  Modified version of \code{s} originally from package \code{mgcv}, Simon Wood (2014).
#' @author Alastair Rushworth
#' @seealso \code{\link[=smnet]{smnet}}, \code{\link{plot.smnet}}, \code{\link{predict.smnet}}
#' @importFrom stats terms
#' @importFrom stats reformulate
#' @export
#' 

m <- function (...,  k = -1, cyclic = F) {
  vars   <- as.list(substitute(list(...)))[-1]
  d      <- length(vars)
  term   <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  
  
  if (d > 1) for (i in 2:d) term[i] <- deparse(vars[[i]], backtick = TRUE, width.cutoff = 500)
  for (i in 1:d)            term[i] <- attr(terms(reformulate(term[i])), "term.labels")
  
  # some warning messages 
  if (all.equal(round(k, 0), k) != TRUE) stop("argument k of m() should be integer")
  if (length(unique(term)) != d)         stop("Repeated variables as arguments of a smooth are not permitted")
  
  
  full.call <- paste("m(", term[1], sep = "")
  if (d > 1) for (i in 2:d) full.call <- paste(full.call, ",", term[i], sep = "")
  label <- paste(full.call, ")", sep = "")
  list(term = term, bs.dim = k,  cyclic = cyclic)
}