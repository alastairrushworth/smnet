\name{plot.smnet}
\alias{plot.smnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot a Stream Network Model}
\description{Plot linear, univariate and bivariate smooth effects and network smooth terms that resulting from a call to \code{\link{smnet}}.  }
\usage{
  \method{plot}{smnet}(x, type = "covariates", se = FALSE, res = FALSE, weight = NULL, 
  sites = FALSE, sites.col = NULL, sites.cex = 1, network.col = NULL,
  shadow = 0, key = TRUE, legend.text = NULL, legend.range = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \bold{smnet}
  }
  \item{type}{Character string identifying the type of plot to be produced.  The default, "\code{covariates}", produces plots of all linear and smooth components (the latter corresponding to each appearance of \code{\link{m}} in the model formula).  "\code{full}" plots the stream network fitted mean using the full set of spatial points contained in the associated \code{SSN} object.   "\code{segments}" plots the stream network fitted mean using a set of connnected line segments to represent the spatial network, this can be faster for large networks than "\code{full}". 
  }
  \item{se}{Logical.  When TRUE and \code{type = }"\code{covariates}", standard errors are shown on plots of linear and smooth terms.  When \code{type = } "\code{segment}" or "\code{full}" spatial standard errors are plotted.}
  \item{res}{Logical.  When TRUE, partial residuals are shown on plots of linear and smooth component.  Ignored when  \code{type = } "\code{full}" or "\code{segments}"}
  \item{weight}{Positive real number that scales stream segment widths (as determined using Shreve order) to indicate relative size of stream segments.  Ignored when \code{type = }"\code{covariate}".  Currently only "\code{autoShreve}" is supported, defaults to \code{NULL} in which all streams segments are plotted with lines with identical widths.
  }
  \item{sites}{Logical indicating whether locations of observation sites should be added to spatial plots.  Ignored when \code{type = }"\code{covariate}" and defaults to \code{FALSE}. 
  }
  \item{sites.col}{Single colour to plot observation locations.  If not provided, points will be coloured according to the default legend and average observation at each location.}
  \item{sites.cex}{Expansion factor for the size of plotted observation points}
    \item{network.col}{Single colour to represent all stream segments. By default, network is coloured according to fitted values from model. Ignored when \code{type = }"\code{covariate}".}
  \item{shadow}{Positive scalar that adds a black outline to spatial stream segment plots, useful if the colour scale has many light colours.  Ignored when \code{type = }"\code{covariate}" and defaults to 0 (no shadow is drawn). 
  }
   \item{key}{Logical.  Plots a colour legend for network plots.  Ignored when \code{type = "covariates"}.
  }
    \item{legend.text}{Character annotation to add to color scale.  Ignored if \code{key = FALSE} or \code{type = "covariates"}.
  }
  \item{legend.range}{Range of values represented by the color scale.  Ignored if \code{key = FALSE} or \code{type = "covariates"}.}
    \item{...}{Other arguments passed to plot}
}

\examples{

# Set up an SSN object - this part taken 
# from the SSN:::SimulateOnSSN help file
example_network_4<- createSSN(
                     n          = 50,
                     obsDesign  = binomialDesign(200), 
                     predDesign = binomialDesign(50),
                     importToR  = TRUE, 
                     path = paste(tempdir(),"/example_network_4",sep = ""),
                     treeFunction = iterativeTreeLayout
                     )

# plot the simulated network structure with prediction locations
# plot(example_network_4, bty = "n", xlab = "x-coord", ylab = "y-coord")

## create distance matrices, including between predicted and observed
createDistMat(example_network_4, "preds", o.write=TRUE, amongpred = TRUE)

## extract the observed and predicted data frames
observed_data            <- getSSNdata.frame(example_network_4, "Obs")
prediction_data          <- getSSNdata.frame(example_network_4, "preds")

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
               ssn.object      = example_network_4, 
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
                  paste(tempdir(), "/example_network_4", sep = ""), 
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


\author{Alastair Rushworth}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  \seealso{
   \code{\link{predict.smnet}}, \code{\link{summary.smnet}}
  }


