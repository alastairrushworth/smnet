#' @importFrom graphics lcm

plot_segments <- function(
  x, weight, netID, se, sites, sites.col, sites.cex, shadow, legend.text,
  legend.range, network.col, key, ...){    
  
  # break up beta_hat
  getcol    <- function(M) ifelse(ncol(M) == "NULL", 1, ncol(M)) 
  cov.dims  <- lapply(x$internals$X.list, getcol)    
  inds      <- unlist(cov.dims)
  ord       <- order(x$internals$ord)
  cum.inds  <- cumsum(inds)
  n.cov     <- length(x$internals$X.list)
  ests      <- vector("list", length = n.cov)
  ests[[1]] <- x$internals$beta_hat[1]
  for(i in 1:(n.cov-1)) ests[[i+1]] <- x$internals$beta_hat[(cum.inds[i]+1):(cum.inds[i+1])]
  spatial_comp <- ests[[n.cov]] + ests[[1]]
  
  if(se){
    l.covs            <- lapply(x$internals$X.list, ncol)
    component.size    <- l.covs[[length(l.covs)]]
    spatial.n         <- nrow(x$internals$adjacency[[1]])
    empty.X.list      <- lapply(l.covs, make_sparse, nrow = spatial.n)
    empty.X.list[[length(l.covs)]] <- diag.spam(spatial.n)
    new.X             <- Reduce("cbind", empty.X.list)
    left1             <- forwardsolve.spam(x$internals$U, t(new.X), transpose = T)
    left2             <- backsolve.spam(x$internals$U, left1, transpose = T)
    vec               <- x$internals$X.spam %*% left2
    spatial_comp      <- sqrt(colSums(vec*vec)*x$internals$sigma.sq)   
  }
  
  fitted_segment <- spatial_comp[ord]
  # SET UP THE COLOUR SCALE FOR PLOTTING
  nlevels      <- 20
  if(is.null(legend.range)){
    zlim         <- range(fitted_segment, finite = TRUE)
  } else {
    zlim         <- legend.range
  }
  brks         <- pretty(zlim, n = 10)
  nclasses     <- length(brks) - 1 
  nnodes       <- length(fitted_segment)
  colorPalette <- rev(heat.colors(nclasses))
  mar.orig     <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w            <- (3 + mar.orig[2L]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  mar     <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  
  # SET UP THE PLOTTING REGION
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(brks), xaxs = "i", yaxs = "i")

  if(key){
    # PLOT THE LEGEND
    rect(0, brks[-length(brks)], 1, brks[-1L], col = colorPalette)
    axis(4)
    mtext(side = 4, legend.text, line = 2)
    box()
  }
  
  # GET WEIGHT IF REQUESTED
  if(!is.null(weight)) shreve  <- adjacency_to_shreve(adjacency = x$internals$adjacency)[ord]
  default.parameters <- list(xlim = range(x$ssn.object@bbox[1,]), ylim = range(x$ssn.object@bbox[2,]))
  updated.parameters <- modifyList(default.parameters, list(...))
  
  # PLOT THE FULL NETWORK
  mar      <- mar.orig
  mar[4L]  <- 1
  par(mar = mar)
  plot.new()
  do.call(plot.window, updated.parameters)
  data      <- x$ssn.object@data
  netIDinds <- which(data$netID == netID)
  data_rid  <- data$rid
  

  
  # create quantile breaks
  # this is the order of the rids in the adjacency matrix
  lower.breaks <- c(min(fitted_segment, na.rm = T) - .0001, brks)
  upper.breaks <- c(brks, max(fitted_segment, na.rm = T))
  if(!shadow == 0){
    # lay down black background shadow
    for(k in 1:nclasses) {
      for(i in 1:length(netIDinds)){
        z <- netIDinds[i]
        lwd_i <- ifelse(!is.null(weight), weight*(log(shreve[i]) + 1), 1)
        for(j in 1:length(x$ssn.object@lines[[z]])){
          if(fitted_segment[i] > lower.breaks[k] & fitted_segment[i] <= upper.breaks[k]){
            coords <- x$ssn.object@lines[[z]]@Lines[[j]]@coords
            lines(coords[c(1, nrow(coords)), ], col = 1, lwd = lwd_i+shadow)
          }
        }
      }
    }
  }
  # plot the colured lines
  for(k in 1:nclasses) {
    for(i in 1:length(netIDinds)){
      z <- netIDinds[i]
      lwd_i <- ifelse(!is.null(weight), weight*(log(shreve[i]) + 1), 1)
      for(j in 1:length(x$ssn.object@lines[[z]])){
        if(fitted_segment[i] > lower.breaks[k] & fitted_segment[i] <= upper.breaks[k]){
          coords <- x$ssn.object@lines[[z]]@Lines[[j]]@coords
          coli   <- ifelse(is.null(network.col), colorPalette[k], network.col)
          lines(coords[c(1, nrow(coords)), ], col = coli, lwd = lwd_i)
        }
      }
    }
  }
  if(sites){
    data_obs      <- getSSNdata.frame(x$ssn.object)
    data_obs      <- data_obs[data_obs$netID == netID, ]
    locs_char     <- paste(data_obs$NEAR_X, data_obs$NEAR_Y, sep = "_")
    mn_resp       <- tapply(x$internals$response, locs_char, mean)
    locs_char     <- names(mn_resp)
    pts_obs       <- Reduce("rbind", lapply(strsplit(locs_char, "[_]"), as.numeric))
    if(is.null(sites.col)){
      col.ints  <- cut(mn_resp, breaks = brks, labels = FALSE, right = FALSE)
      sites.col <- colorPalette[col.ints]
    }
    points(pts_obs[,1], pts_obs[,2], pch = 20, col = 1, cex = 1.5*sites.cex)
    points(pts_obs[,1], pts_obs[,2], pch = 20, col = sites.col, cex = sites.cex)
  }
  do.call(title, updated.parameters)
  if(is.null(updated.parameters$xaxt)){
    xlimits <- updated.parameters$xlim
    xticks  <- format(seq(xlimits[1], xlimits[2], length.out = 5), digits = 2, zero.print = T)
    Axis(x = xlimits, at = xticks, side = 1, labels = xticks)
  }
  if(is.null(updated.parameters$yaxt)){
    ylimits <- updated.parameters$ylim
    yticks  <- format(seq(ylimits[1], ylimits[2], length.out = 5), digits = 2, zero.print = T)
    Axis(x = ylimits, at = yticks, side = 2, labels = yticks)
  }
}

