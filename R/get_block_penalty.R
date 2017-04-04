
# utility to calculate the block diagonal penalty matrices for spline penalisation

get_block_penalty<-function(P, blockzero, i){        
  blockzero[[i]] <- P
  b.diag         <- function(L){
    diag.spam <- L[[1]]
    for(j in 2:length(L)){
      diag.spam <- bdiag.spam(diag.spam, L[[j]])
    }
    diag.spam
  }
  b.diag(blockzero)
}	
