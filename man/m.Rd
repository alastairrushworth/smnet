\name{m}
\alias{m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specify Smooth Terms in Formulae
}
\description{
Function used to set up univariate or bivariate smooth terms based on P-splines, for use within a call to \code{\link{smnet}}.   
}
\usage{
m(..., k = -1, cyclic = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
one or more variables for creating P-spline smooths
}
  \item{k}{
integer defining the number of uniformly spaced B-spline basis functions for the smooth, default is 10.  For 2d (and higher) smooths, this is the marginal basis size.
}
 \item{cyclic}{
logical vector indicating whether the smooth should be cyclic.  Based on the harmonic smoother of Eiler and Marx (2004)
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{term}{character vector of the names of the variables involved in the smooth to be set up}
\item{bs.dim}{number of B-spline basis functions to be used in the smooth}
}
\references{
Modified version of \code{s} originally from package \code{mgcv}, Simon Wood (2014).
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
example_network_2<- createSSN(
                     n          = 50,
                     obsDesign  = binomialDesign(200), 
                     predDesign = binomialDesign(50),
                     importToR  = TRUE, 
                     path = paste(tempdir(),"/example_network_2",sep = ""),
                     treeFunction = iterativeTreeLayout
                     )

# plot the simulated network structure with prediction locations
# plot(example_network_2, bty = "n", xlab = "x-coord", ylab = "y-coord")

## create distance matrices, including between predicted and observed
createDistMat(example_network_2, "preds", o.write=TRUE, amongpred = TRUE)

## extract the observed and predicted data frames
observed_data            <- getSSNdata.frame(example_network_2, "Obs")
prediction_data          <- getSSNdata.frame(example_network_2, "preds")

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
               ssn.object      = example_network_2, 
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
                  paste(tempdir(), "/example_network_2", sep = ""), 
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

