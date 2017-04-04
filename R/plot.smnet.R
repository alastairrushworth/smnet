plot.smnet<- function(x, 
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
  
  # ### --------------------------------------------
  # ### PLOT NETWORK ON A GOOGLE MAP  - coming soon
  # ### --------------------------------------------
  # if(type == "gmaps"){
  # googleOptions = NULL,
  #   zoom        <- googleOptions$zoom
  #   SCALE       <- googleOptions$SCALE
  #   maptype     <- googleOptions$maptype
  #   destfile    <- ifelse(is.null(googleOptions$destfile),  "map.png", googleOptions$destfile)
  #   if(!x$internals$net) stop("No spatial network component to plot in x")
  #   plot_gmaps(x, weight = weight, netID = netID, se = se, sites = sites, 
  #              shadow = shadow, legend.text = legend.text, key = key,
  #              maptype = maptype, zoom = zoom, SCALE = SCALE,
  #              sites.col = sites.col, sites.cex = sites.cex, network.col = network.col, 
  #              destfile = destfile, ...)
  # }
  
  ### ----------------------------------------
  ### PLOT SMOOTH AND LINEAR COVARIATES 
  ### ----------------------------------------
  if(type == "covariates"){
    plot_cov(x = x$internals, res = res, ...)
  }
}
 


