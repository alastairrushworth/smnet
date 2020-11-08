#' @importFrom spam kronecker.spam

construct_bspline_penalty <- function(variables_length, dimensions, term.number, blockzero, cyclic){
  # if only a single variable, then build spline basis
  if(variables_length == 1){
    P              <- get_penalty(dimensions, cyclic = cyclic)
    P_wee_list     <- list(P = P)
    P_list         <- list(P = get_block_penalty(P, blockzero = blockzero, i = term.number))
  } else if(variables_length >  1){
    # if more than one variable put each resulting penalty in a list, then Reduce using boxproduct
    if(length(dimensions) == 1){
      dim.each <- rep(dimensions, variables_length)
    } else if(length(dimensions) == variables_length){
      dim.each <- dimensions
    } else stop("the vector `k' is the wrong length for one of the m() terms")
    if(length(cyclic) == 1){
      cyclic.each <- rep(cyclic, variables_length)
    } else if(length(cyclic) == variables_length){
      cyclic.each <- cyclic
    } else stop("provided logical vector `cyclic' is the wrong length for one of the m() terms")
    marginal_penalty_list <- vector("list", length = variables_length)
    for(i in 1:variables_length){
      marginal_penalty_list[[i]] <- get_penalty(dim.each[i], cyclic = cyclic.each[i])
    }
    # need to find way to construct all combinations of kroneckers with identity matrix
    kronecker_list   <- vector("list", length = variables_length)
    marginal_id_list <- vector("list", length = variables_length)
    for(i in 1:variables_length) marginal_id_list[[i]] <- diag.spam(1, dim.each[i]) 
    for(i in 1:variables_length){
      int_list            <- marginal_id_list
      int_list[[i]]       <- marginal_penalty_list[[i]]
      kronecker_list[[i]] <- Reduce("kronecker.spam", int_list)
    }
    P_list     <- lapply(kronecker_list, get_block_penalty, blockzero = blockzero, i = term.number)
    P_wee_list <- kronecker_list
  }
  list(P_list = P_list, P_wee_list = P_wee_list)
}
