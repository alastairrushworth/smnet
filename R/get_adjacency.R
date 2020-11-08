#' Construct a Stream Network Adjacency Matrix
#' 
#' @description Builds a sparse adjacency matrix from a user specified 
#' \code{SSN} data directory, by extracting and processing the 
#' \code{binaryID.db} table.  The resulting output of this function is 
#' required input for fitting spatial additive network models to \code{SSN} 
#' objects using the main \code{\link{smnet}} function.
#' 
#' @param ssn_directory Required character string indicating the path to the 
#' location of the .ssn directory which contains the binaryID.db table.
#' @param netID Logical.  Integer specifying the particular stream network 
#' of interest within the \code{SSN} object.  Defaults to 1.
#' 
#' @return List object with components
#' \itemize{
#' \item{\code{adjacency}: Sparse adjacency matrix of class \code{spam} with row and 
#' column dimension equal to the number of stream segments.  If the i^th column 
#' has non-zero elements \eqn{j_1}{j_1} and \eqn{j_2}{j_2} then this indicates 
#' that \eqn{j_1}{j_1} and \eqn{j_2}{j_2} are direct upstream neighbours of i.
#'   If the \eqn{i^\textrm{th}}{i^th} column has sum 1, then this indicates that 
#'   \eqn{i}{i} has only one upstream neighbour, and therefore no confluence lies
#'    between them; by default the spatial penalties treat these differently.}
#' \item{\code{bid}: Character vector of binary identifiers for each stream segment, 
#' used only for automatic calculation of Shreve's stream order within 
#' \code{smnet}}
#' }
#' 
#' @author Alastair Rushworth
#' @seealso   \code{\link[=smnet]{smnet}}, \code{\link{plot.smnet}}, \code{\link{predict.smnet}}
#' @export
#' @importFrom spam spam

get_adjacency <- function(ssn_directory, netID = 1){
  binaryIDs  <- get_binaryIDs(ssn_directory, net = netID)
  bid        <- binaryIDs[, 2]
  rid        <- binaryIDs[, 1]
  nch        <- nchar(bid)
  bid.list   <- split(bid, nch)
  rid.list   <- split(rid, nch)
  n.segments <- nrow(binaryIDs)
  xy         <- matrix(ncol = 2, nrow = (n.segments - 1))
  counter    <- 0
  for(j in 1:(length(bid.list) - 1)){
    bid.dn  <- bid.list[[j]]
    bid.up  <- bid.list[[j + 1]]
    rid.dn  <- rid.list[[j]]
    rid.up  <- rid.list[[j + 1]]
    bid.up.sub.vec <- bid.up
    rid.up.sub.vec <- rid.up
    for(i in 1:length(bid.dn)){
      current.dn.bid  <- bid.dn[i]
      current.dn.rid  <- rid.dn[i]
      inner.count     <- 1
      number.upstream <- 0
      n.bid.up        <- length(bid.up.sub.vec)
      crit            <- FALSE
      while(!crit){	
        if(n.bid.up > 0){
          current.up.bid <- bid.up.sub.vec[inner.count]
          connected <- substr(current.up.bid, 1, nchar(current.dn.bid)) == current.dn.bid
          if(connected){
            counter         <- counter + 1
            number.upstream <- number.upstream + 1
            xy[counter, ]   <- c(current.dn.rid, rid.up.sub.vec[inner.count])
            rid.up.sub.vec  <- rid.up.sub.vec[-inner.count]
            bid.up.sub.vec  <- bid.up.sub.vec[-inner.count]
          }
          if(!connected) inner.count <- inner.count + 1
          crit <- (number.upstream == 2) | ((number.upstream + inner.count - 1) == n.bid.up)
        }
        if(n.bid.up == 0) crit <- TRUE
      }
    }		
  }	
  xy[,1]  <- re_map_rid(rid_vector = xy[, 1], all_rid = rid)
  xy[,2]  <- re_map_rid(rid_vector = xy[, 2], all_rid = rid)
  add.one <- min(xy) == 0 # shouldn't be needed if remapping worked properly
  list(
    adjacency = spam(
      list(j = (xy[, 1] + add.one), 
           i = (xy[, 2] + add.one), 
           rep(1, nrow(xy))), 
      nrow = n.segments, 
      ncol = n.segments
    ), 
    rid_bid = cbind(rid, bid)
  )
}	

  