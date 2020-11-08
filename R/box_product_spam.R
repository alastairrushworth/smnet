
not_sparse_box_product<-function(A, B){
  unit.vec.B <- matrix(1, nrow = 1, ncol = ncol(B))
  unit.vec.A <- matrix(1, nrow = 1, ncol = ncol(A))
  (A %x% unit.vec.B) * (unit.vec.A %x% B)}
