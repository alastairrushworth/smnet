context('Testing lsn')

mf_location <- "/Users/alastairrushworth/Documents/git_repositories/smnet/tests/testthat/testdata/lsn.ssn"

mf_location <- paste0(normalizePath('testdata'), "/lsn.ssn")
mf          <- suppressWarnings(SSN::importSSN(mf_location, predpts = "preds", o.write = TRUE))
adjacency   <- get_adjacency(mf_location, netID = 1)


# test that maxit works properly
test_that("Linear model with network", {
  # simple gam
  mod_smn     <- 
    smnet(
      formula = ELEV ~ CUMDRAINAG + network(adjacency = adjacency, weight = "areaPI"),
      data.object = mf, netID = 1, control = list(verbose = FALSE))
  expect_equal(class(mod_smn), 'smnet')
  expect_equal(names(summary(mod_smn, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  plot(mod_smn, type = "covariates")
  plot(mod_smn, type = "segments", weight = 0.2, shadow = 2)
  plot(mod_smn, type = "full", weight = 0.2, shadow = 2)
  expect_equal(
    names(predict(mod_smn)),
    c('predictions', 'predictions.se')
  )
  expect_equal(
    names(predict(mod_smn, newdata = getSSNdata.frame("lsn.ssn", "Obs"))),
    c('predictions', 'predictions.se')
  )
})






