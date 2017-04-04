\name{get_adjacency}
\alias{get_adjacency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct an Adjacency Matrix}
\description{Builds a sparse adjacency matrix from a user specified \code{SSN} data directory, by extracting and processing the binaryID.db table.  The resulting output of this function is required input for fitting spatial additive network models to \code{SSN} objects using the main \code{\link{smnet}} function.}
\usage{
get_adjacency(ssn_directory, netID)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ssn_directory}{Required character string indicating the path to the location of the .ssn directory which contains the binaryID.db table}
   \item{netID}{Integer specifying the particular stream network of interest within the \code{SSN} object.  Defaults to 1.}
}

\value{List with two elements.  
\item{adjacency}{Sparse adjacency matrix of class "\code{spam}" with row and column dimension equal to the number of stream segments.  If the i^th column has non-zero elements \eqn{j_1}{j_1} and \eqn{j_2}{j_2} then this indicates that \eqn{j_1}{j_1} and \eqn{j_2}{j_2} are direct upstream neighbours of i.  If the \eqn{i^\textrm{th}}{i^th} column has sum 1, then this indicates that \eqn{i}{i} has only one upstream neighbour, and therefore no confluence lies between them; by default the spatial penalties treat these differently.}

\item{bid}{Character vector of binary identifiers for each stream segment, used only for automatic calculation of Shreve's stream order within \link{smnet}}
}

\author{Alastair Rushworth}



\seealso{
\code{\link{smnet}}
}
\examples{
# Set up an SSN object - this part taken 
# from the SSN:::SimulateOnSSN help file
example_network_1<- createSSN(
                     n          = 50,
                     obsDesign  = binomialDesign(200), 
                     predDesign = binomialDesign(50),
                     importToR  = TRUE, 
                     path = paste(tempdir(),"/example_network_1",sep = ""),
                     treeFunction = iterativeTreeLayout
                     )

# plot the simulated network structure with prediction locations
# plot(example_network, bty = "n", xlab = "x-coord", ylab = "y-coord")

## create distance matrices, including between predicted and observed
createDistMat(example_network_1, "preds", o.write=TRUE, amongpred = TRUE)

## extract the observed and predicted data frames
observed_data            <- getSSNdata.frame(example_network_1, "Obs")
prediction_data          <- getSSNdata.frame(example_network_1, "preds")

## associate continuous covariates with the observation locations
#  data generated from a normal distribution
obs                      <- rnorm(200)
observed_data[,"X"]      <- obs
observed_data[,"X2"]     <- obs^2

## associate continuous covariates with the prediction locations
#  data generated from a normal distribution
pred                     <- rnorm(50) 
prediction_data[,"X"]    <- pred
prediction_data[,"X2"]   <- pred^2

## simulate some Gaussian data that follows a 'tail-up' spatial process
sims <- SimulateOnSSN(
               ssn.object      = example_network_1, 
               ObsSimDF        = observed_data, 
               PredSimDF       = prediction_data,  
               PredID          = "preds",  
               formula         = ~ 1 + X,
               coefficients    = c(1, 10),
               CorModels       = c("Exponential.tailup"), 
               use.nugget      = TRUE,
               CorParms        = c(10, 5, 0.1),
               addfunccol      = "addfunccol")$ssn.object


## extract the observed and predicted data frames, now with simulated values
sim1DFpred         <- getSSNdata.frame(sims, "preds")
sim1preds          <- sim1DFpred[,"Sim_Values"]
sim1DFpred[,"Sim_Values"] <- NA
sims               <- putSSNdata.frame(sim1DFpred, sims, "preds")

# create the adjacency matrix for use with smnet
adjacency    <- get_adjacency(
                  paste(tempdir(), "/example_network_1", sep = ""), 
                  net = 1
                  )

# not run - plot the adjacency matrix
# display(adjacency[[1]])

# sometimes it is useful to see which varaibles are valid network weights 
# in the data contained within the SSN object
show_weights(sims, adjacency)

# fit a penalised spatial model to the stream network data
# Sim_Values are quadratic in the X covariate.  To highlight 
# the fitting of smooth terms, this is treated as non-linear 
# and unknown using m().
mod_smn       <- smnet(formula = Sim_Values ~ m(X) + m(X2) + 
network(adjacency = adjacency, weight = "shreve"), 
                       data.object = sims, netID = 1)

# not run - plot different summaries of the model
plot(mod_smn, type = "covariates")
plot(mod_smn, type = "segments", weight = 4, shadow = 2)
plot(mod_smn, type = "full", weight = 4, shadow = 2)

# obtain predictions at the prediction locations and plot 
# against true values
preds <- predict(mod_smn, newdata = getSSNdata.frame(sims, "preds"))
plot(preds$predictions, sim1preds)

# obtain summary of the fitted model
summary(mod_smn)
}





% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{network}
\keyword{sparse}
\keyword{P-spline}