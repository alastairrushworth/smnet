\name{summary.smnet}
\alias{summary.smnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarise Stream Network Model
}
\description{
Generate summaries of linear and smooth components of an \code{smnet} object.
}
\usage{
\method{summary}{smnet}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{smnet}.
}
\item{...}{
other arguments passed to summary
}
}

\value{List object with components
\item{1}{\code{linear.terms}: the linear components of the fitted model}
\item{2}{\code{smooth.terms}: the values of the smoothing parameters on the log scale, and the partial degrees of freedom associated with each smooth component.  Note: Network components always have two smoothing parameters, where the second is a (usually small) ridge parameter.}
}


\author{
Alastair Rushworth
}

\seealso{
\code{\link{smnet}}
}
\examples{

# Set up an SSN object - this part taken 
# from the SSN:::SimulateOnSSN help file
example_network_9<- createSSN(
                     n          = 50,
                     obsDesign  = binomialDesign(200), 
                     predDesign = binomialDesign(50),
                     importToR  = TRUE, 
                     path = paste(tempdir(),"/example_network_9",sep = ""),
                     treeFunction = iterativeTreeLayout
                     )

# plot the simulated network structure with prediction locations
# plot(example_network_9, bty = "n", xlab = "x-coord", ylab = "y-coord")

## create distance matrices, including between predicted and observed
createDistMat(example_network_9, "preds", o.write=TRUE, amongpred = TRUE)

## extract the observed and predicted data frames
observed_data            <- getSSNdata.frame(example_network_9, "Obs")
prediction_data          <- getSSNdata.frame(example_network_9, "preds")

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
               ssn.object      = example_network_9, 
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
                  paste(tempdir(), "/example_network_9", sep = ""), 
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

