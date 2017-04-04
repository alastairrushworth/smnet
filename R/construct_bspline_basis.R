construct_bspline_basis <- function(variables, dimensions, cyclic, range.variables = variables){
  # if only a single variable, then build spline basis
  variables.mat       <- as.matrix(variables)
  range.variables.mat <- as.matrix(range.variables)
  n.varb              <- ncol(variables.mat)
  if(n.varb == 1){
    if(cyclic){
      bspline_basis <- cyclic_b_spline_basis(x = variables.mat[, 1], nseg = (dimensions - 3), 
                                             range.variables = range.variables.mat)
    } else {
      bspline_basis <- b_spline_basis(x = variables.mat[, 1], nseg = (dimensions - 3), 
                                      range.variables = range.variables.mat) 
    }
    bspline_basis <- make_spam(bspline_basis)
  } 
  # if there are several variables for a tensor product smooth, construct model matrix
  if(n.varb > 1){
    # if only one dimension is provided, then recycle it of if the dimension is the wrong length
    if(length(dimensions) == 1){
      dim.each <- rep(dimensions, n.varb)
    } else if(length(dimensions) == n.varb){
      dim.each <- dimensions
    } else stop("the vector `k' is the wrong length for one of the m() terms")
    # if only one cyclic is provided, then recycle it of if the cyclic is the wrong length
    if(length(cyclic) == 1){
      cyclic.each <- rep(cyclic, n.varb)
    } else if(length(cyclic) == n.varb){
      cyclic.each <- cyclic
    } else stop("provided logical vector `cyclic' is the wrong length for one of the m() terms")
    basis_list <- vector("list", length = n.varb)
    # loop through variables and create marginal basis matrices
    for(i in 1:n.varb){
      if(cyclic.each[i]){
        basis_list[[i]] <- cyclic_b_spline_basis(x = variables.mat[, i], nseg = dim.each[i] - 3, 
                                                 range.variables = range.variables.mat[, i])
      } else {
        basis_list[[i]] <- b_spline_basis(x = variables.mat[, i], nseg = dim.each[i] - 3, 
                                          range.variables = range.variables.mat[, i])
      }
    }
    # calculated the row-wise tensor product of the marginal basis matrices
    bspline_basis <- make_spam(Reduce("not_sparse_box_product", basis_list))
  }
  bspline_basis
}

