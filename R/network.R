#' @export
network <- function(adjacency = NULL,  weight = "autoShreve", fixed.df = NULL){
  ret <- list(adjacency = adjacency,  weight = weight, fixed.df = fixed.df, by = "NA")
  class(ret) <- "network.spec"
  ret
}
