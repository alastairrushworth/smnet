\name{show_weights}
\alias{show_weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Search for and Validate Weights in an SSN Object}

\description{Explore SSN objects for valid stream weights for use in fitting stream network models.}
\usage{
show_weights(SSNobject, adjacency, netID = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SSNobject}{\code{SpatialStreamNetwork} containing data to be searched for valid network weights}
  \item{adjacency}{adjacency object corresponding to \code{SSNobject}, resulting from a call to \link{get_adjacency}}
  \item{netID}{Positive integer indicating the network number to investigate, if multiple networks are contained in\code{SSNobject}.  Default is 1.}
}

\value{Prints the names of valid weight variables to the console.}

\author{Alastair Rushworth}


\seealso{
\link{get_adjacency}, \link{network}
}
\examples{

# Set up an SSN object - this part taken 
# from the SSN:::SimulateOnSSN help file
example_network_6<- createSSN(
                     n          = 50,
                     obsDesign  = binomialDesign(200), 
                     predDesign = binomialDesign(50),
                     importToR  = TRUE, 
                     path = paste(tempdir(),"/example_network_6",sep = ""),
                     treeFunction = iterativeTreeLayout
                     )

# plot the simulated network structure with prediction locations
# plot(example_network_6, bty = "n", xlab = "x-coord", ylab = "y-coord")

## create distance matrices, including between predicted and observed
createDistMat(example_network_6, "preds", o.write=TRUE, amongpred = TRUE)

## extract the observed and predicted data frames
observed_data            <- getSSNdata.frame(example_network_6, "Obs")
prediction_data          <- getSSNdata.frame(example_network_6, "preds")

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
               ssn.object      = example_network_6, 
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
                  paste(tempdir(), "/example_network_6", sep = ""), 
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

