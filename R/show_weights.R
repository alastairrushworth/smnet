#' Search for and Validate Weights Columns in an SSN Object
#' 
#' @description Explore SSN objects for valid stream weights for use in fitting stream network models.
#' 
#' @param SSNobject \code{SpatialStreamNetwork} containing data to be searched for valid network weights
#' for the smooth, default is 10.  For 2d (and higher) smooths, this is the marginal 
#' basis size.
#' @param adjacency adjacency object corresponding to \code{SSNobject}, resulting from a call to \code{get_adjacency()}.
#' @param netID Positive integer indicating the network number to investigate, if multiple networks are contained in \code{SSNobject}.  Default is 1.
#'
#'  @return Prints the names of valid weight variables to the console.
#' 
#' @author Alastair Rushworth
#' @seealso   \code{\link[=smnet]{smnet}}, \code{\link{plot.smnet}}, \code{\link{predict.smnet}}
#' @export
#' 
show_weights <- function(SSNobject, adjacency, netID = 1){
  # extract the network structure data 
  rid_data        <- SSNobject@data
  # arrange the rid's in ascending order
  ord             <- order(rid_data$rid)
  rid_data        <- rid_data[ord,]
  # extract only the rid's correspoding to the selected netID
  rid_data        <- rid_data[rid_data$netID == netID,]
  # extract the variable names in the data object
  variable_names  <- names(rid_data)
  # vector to store the derived weight types in - ie, whether valid or not
  weight_type     <- vector("character", length = length(variable_names))
  # now cycle through columns in turn and work out which is a valid weight
  for(i in 1:length(variable_names)){
    # this is the info contained in the i^th column
    current_var <- as.vector(rid_data[, i])
    # possibly the data could be non-numeric, in which case skip this part
    if((class(current_var) %in% c("integer", "numeric")) & (!all(current_var == 0))){
      # if the data are numeric, then try to check the weight type
      weight_type[i] <- check_weight(wgt = current_var, adj = adjacency) 
    } else {
      # if non-numeric, just assign unrecognised
      weight_type[i] <- "unrecognised"
    }
  }
  weight.poss     <- variable_names[which(!weight_type == "unrecognised")]
  weight.poss     <- paste(paste(weight.poss, collapse = ", "), ".", sep = "")
  weight.poss     <- paste("\nThe following recognised weights were found: \n", 
                           "-------------------------------------------- \n", 
                           weight.poss, sep = "")
  cat(weight.poss)
}