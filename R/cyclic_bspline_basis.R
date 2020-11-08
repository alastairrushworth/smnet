#' @importFrom splines splineDesign
# Construct cyclic sparse B-spline basis
cyclic_b_spline_basis  	<-	function(x, nseg = 10, deg = 3, range.variables){
  nseg            <- nseg + 3
  x_min           <- min(range.variables)
  x_max           <- max(range.variables)
  interval        <- (x_max - x_min) /nseg
  knot_seq        <- seq(x_min - deg * interval, x_max + deg * interval, by = interval)
  sparseBasis     <- splineDesign(knots = knot_seq, x = x, ord = deg + 1, outer.ok = T, sparse = T)
  lastdegcols     <- (ncol(sparseBasis) - deg + 1):ncol(sparseBasis)
  make_spam(cbind((sparseBasis[, 1:deg] + sparseBasis[, lastdegcols]), sparseBasis[, ((deg + 1):(ncol(sparseBasis) - deg)) ]))
}
