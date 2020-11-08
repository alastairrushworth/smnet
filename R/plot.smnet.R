#' Visualise a Fitted Stream Network Model
#' 
#' @description Plot linear, univariate and bivariate smooth effects and network 
#' smooth terms that resulting from a call to \code{smnet()}.
#' 
#' @param x An object of class \bold{smnet}.
#' @param type Character string identifying the type of plot to be produced.  
#' The default, \code{type = 'covariates'}, produces plots of all linear and smooth 
#' components (the latter corresponding to each appearance of \code{m()} 
#' in the model formula).  \code{type = 'full'} plots the stream network fitted mean 
#' using the full set of spatial points contained in the associated \code{SSN} 
#' object.   \code{type = 'segments'} plots the stream network fitted mean using a set 
#' of connected line segments to represent the spatial network, this can be 
#' faster for large networks than \code{type = 'full'}.
#' @param se Logical.  When \code{TRUE} and \code{type = 'covariates'}, standard errors 
#' are shown on plots of linear and smooth terms.  When \code{type = 'segment'} or 
#' \code{type = 'full'} spatial standard errors are plotted.
#' @param res Logical.  When \code{TRUE}, partial residuals are shown on plots of linear 
#' and smooth component.  Ignored when  \code{type = 'full'} or \code{type = 'segments'}.
#' @param weight Positive real number that scales stream segment widths (as determined using 
#' Shreve order) to indicate relative size of stream segments.  Ignored when 
#' \code{type = 'covariate'}.  Currently only \code{'autoShreve'} is supported, defaults to 
#' \code{NULL} in which all streams segments are plotted with lines with identical widths.
#' @param sites Logical indicating whether locations of observation sites should be added to 
#' spatial plots.  Ignored when \code{type = 'covariate'} and defaults to \code{FALSE}.
#' @param sites.col Single colour to plot observation locations.  If not provided, points will 
#' be coloured according to the default legend and average observation at each location.
#' @param sites.cex Expansion factor for the size of plotted observation points.
#' @param network.col Single colour to represent all stream segments. By default, network is
#'  coloured according to fitted values from model. Ignored when \code{type = 'covariate'}.
#' @param shadow Positive scalar that adds a black outline to spatial stream segment plots, 
#' useful if the colour scale has many light colours.  Ignored when \code{type = 'covariate'} 
#' and defaults to 0 (no shadow is drawn).
#' @param key Logical.  Plots a colour legend for network plots.  Ignored when 
#' \code{type = 'covariates'}.
#' @param legend.text Character annotation to add to color scale.  Ignored if 
#' \code{key = FALSE} or \code{type = 'covariates'}.
#' @param legend.range Range of values represented by the color scale.  Ignored if 
#' \code{key = FALSE} or \code{type = 'covariates'}.
#' @param ... Other arguments passed to \code{plot()}.
#' 
#' @author Alastair Rushworth
#' @export

plot.smnet<- function(
  x, 
  type   = "covariates", 
  se     = FALSE, 
  res    = FALSE, 
  weight = NULL,
  sites  = FALSE,
  sites.col   = NULL,
  sites.cex   = 1,
  network.col = NULL,
  shadow = 0,
  key    = TRUE,
  legend.text = NULL,
  legend.range = NULL,
  ...)
{
  ### ----------------------------------------
  ### DO SIMPLE CHECK ON THE INPUT OBJECT
  ### ----------------------------------------
  netID    <- x$internals$netID
  if(!class(x) == "smnet") stop("x must be an object of class 'smnet'")
  if((type %in% c("network-segments", "network-full", "network-gmaps")) && (res == T)) warning("Ignoring 'res' argument, since spatial plot requested")
  

  ### -------------------------------------------
  ### PLOT NETWORK USING CONNECTED LINE SEGMENTS
  ### -------------------------------------------
  if(type == "segments"){
    if(!x$internals$net) stop("No spatial network component to plot in x")
    plot_segments(x, weight = weight, netID = netID, se = se, sites = sites, shadow = shadow, 
                  legend.text = legend.text, legend.range = legend.range, key = key, 
                  sites.col = sites.col, sites.cex = sites.cex,
                  network.col = network.col, ...) 
  }
  
  ### --------------------------------------------
  ### PLOT NETWORK USING ALL AVAILABLE POINTS INFO 
  ### --------------------------------------------
  if(type == "full"){
    if(!x$internals$net) stop("No spatial network component to plot in x")
      plot_full(x, weight = weight, netID = netID, se = se, sites = sites, shadow = shadow, 
                legend.text = legend.text, legend.range = legend.range, key = key, 
                sites.col = sites.col, sites.cex = sites.cex,
                network.col = network.col, ...) 
  }
  
  ### ----------------------------------------
  ### PLOT SMOOTH AND LINEAR COVARIATES 
  ### ----------------------------------------
  if(type == "covariates"){
    plot_cov(x = x$internals, res = res, ...)
  }
}
 


