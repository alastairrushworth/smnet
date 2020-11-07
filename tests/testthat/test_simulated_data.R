

context('Testing simulated data: summary and predict functions')
# Set up an SSN object - this part taken 
# from the SSN:::SimulateOnSSN help file
set.seed(2015)
example_network<- 
  SSN::createSSN(
    n          = 100,
    obsDesign  = binomialDesign(200), 
    predDesign = binomialDesign(50),
    importToR  = TRUE, 
    path = paste(tempdir() ,"/example_network", sep = "")
  )

## create distance matrices, including between predicted and observed
SSN::createDistMat(
  example_network, 
  "preds", 
  o.write = TRUE, 
  amongpred = TRUE
)

## extract the observed and predicted data frames
observed_data            <- SSN::getSSNdata.frame(example_network, "Obs")
prediction_data          <- SSN::getSSNdata.frame(example_network, "preds")

## associate continuous covariates with the observation locations
#  data generated from a normal distribution
obs1                      <- rnorm(200)
obs2                      <- rnorm(200)
obs3                      <- rnorm(200)
observed_data[,"X1"]      <- obs1
observed_data[,"X2"]      <- obs2
observed_data[,"X3"]      <- obs3

## associate continuous covariates with the prediction locations
#  data generated from a normal distribution
pred1                     <- rnorm(50) 
pred2                     <- rnorm(50) 
pred3                     <- rnorm(50) 
prediction_data[,"X1"]    <- pred1
prediction_data[,"X2"]    <- pred2
prediction_data[,"X3"]    <- pred3

## simulate some Gaussian data that follows a 'tail-up' spatial process
sims <- SSN::SimulateOnSSN(
  ssn.object      = example_network, 
  ObsSimDF        = observed_data, 
  PredSimDF       = prediction_data,  
  PredID          = "preds",  
  formula         = ~ 1 + X1 + X2 + X3,
  coefficients    = c(1, 2, 3, 4),
  CorModels       = c("Exponential.tailup"), 
  use.nugget      = TRUE,
  CorParms        = c(10, 5, 0.1),
  addfunccol      = "addfunccol")$ssn.object

## extract the observed and predicted data frames, now with simulated values
sim1DFpred         <- SSN::getSSNdata.frame(sims, "preds")
sim1preds          <- sim1DFpred[,"Sim_Values"]
sim1DFpred[,"Sim_Values"] <- NA
sims               <- SSN::putSSNdata.frame(sim1DFpred, sims, "preds")

# create the adjacency matrix for use with smnet
adjacency    <- get_adjacency(
  paste(tempdir(), "/example_network", sep = ""), 
  net = 1
)

test_that("Single variable GAM", {
  # simple gam
  mod_smn       <- 
    smnet(formula = Sim_Values ~ m(X1, k = 15), 
          data.object = sims, 
          netID = 1, 
          control = list(verbose = FALSE))
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
})


test_that("Single variable cyclic GAM", {
  mod_smn       <- 
    smnet(formula = Sim_Values ~ m(X1, k = 15, cyclic = T), 
          data.object = sims, 
          netID = 1, 
          control = list(verbose = FALSE))
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
  # plot(mod_smn, type = "covariates")
})

test_that("Bivariate GAM", {
  mod_smn       <- 
    smnet(formula = Sim_Values ~ m(X1, X2, k = c(10, 20)),
          data.object = sims, 
          netID = 1, 
          control = list(verbose = FALSE))
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
  # plot(mod_smn, type = "covariates")
})

test_that("Bivariate cyclic GAM", {
  mod_smn       <- 
    smnet(formula = Sim_Values ~ m(X1, X2, k = c(15, 10), cyclic = c(T, T)), 
          data.object = sims, 
          netID = 1, 
          control = list(verbose = FALSE))
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
  # plot(mod_smn, type = "covariates")
})

test_that("3D GAM", {
  mod_smn       <- 
    smnet(formula = Sim_Values ~ X3 + m(X1, X2), 
          data.object = sims, 
          netID = 1, 
          control = list(verbose = FALSE))
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
  # plot(mod_smn, type = "covariates")
})

test_that("Anova GAM", {
  mod_smn       <- 
    smnet(formula = Sim_Values ~ m(X1) + m(X2) + m(X1, X2), 
          data.object = sims, 
          netID = 1, 
          control = list(verbose = FALSE))
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
  # plot(mod_smn, type = "covariates")
})

test_that("Maxit 100", {
  # test that maxit works properly
  mod_smn       <- smnet(
    formula = Sim_Values ~ network(adjacency = adjacency, weight = "shreve"), 
    data.object = sims, 
    netID = 1, 
    control = list(maxit = 100, verbose = FALSE)
  )
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
})

test_that("Maxit 10", {
  # test that maxit works properly
  mod_smn       <- smnet(
    formula = Sim_Values ~ network(adjacency = adjacency, weight = "shreve"), 
    data.object = sims, 
    netID = 1, 
    control = list(maxit = 10, verbose = FALSE)
  )
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
})

test_that("tol <0.000001", {
  # test that tol works properly
  mod_smn       <- smnet(
    formula = Sim_Values ~ network(adjacency = adjacency, fixed.df = NULL,weight = "shreve"), 
    data.object = sims, 
    netID = 1, 
    control = list(tol = 10^-10, verbose = FALSE)
  )
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
})

test_that("tol 0.1", {
  # test that maxit works properly
  mod_smn       <- smnet(
    formula = Sim_Values ~ network(adjacency = adjacency, fixed.df = NULL, weight = "shreve"), 
    data.object = sims, 
    netID = 1, 
    control = list(tol = 10^-1, verbose = FALSE)
  )
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
})




test_that("Fixed df: 10", {

  mod_smn       <- smnet(
    formula = Sim_Values ~ network(adjacency = adjacency, fixed.df = 10, weight = "shreve"), 
    data.object = sims, 
    netID = 1, 
    control = list(verbose = FALSE)
  )
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
  x <- summary(mod_smn, verbose = FALSE)$smooth
  expect_equal(round(as.numeric(as.character(x[[3]])), 0), 10)
})

test_that("Fixed df: 5", {
  
  mod_smn       <- smnet(
    formula = Sim_Values ~ network(adjacency = adjacency, fixed.df = 5, weight = "shreve"), 
    data.object = sims, 
    netID = 1, 
    control = list(verbose = FALSE)
  )
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
  x <- summary(mod_smn, verbose = FALSE)$smooth
  expect_equal(round(as.numeric(as.character(x[[3]])), 0), 5)
})

test_that("Fixed df: 1", {
  
  mod_smn       <- smnet(
    formula = Sim_Values ~ network(adjacency = adjacency, fixed.df = 1, weight = "shreve"), 
    data.object = sims, 
    netID = 1, 
    control = list(verbose = FALSE)
  )
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  expect_equal(names(predict(mod_smn)), c("predictions", "predictions.se"))
  x <- summary(mod_smn, verbose = FALSE)$smooth
  expect_equal(round(as.numeric(as.character(x[[3]])), 0), 1)
})
