context('Testing Middlefork: maxit & plots')

mf_location <- paste0(normalizePath('testdata'), "/MiddleFork/MiddleFork.ssn")
mf          <- importSSN(mf_location, predpts = NULL, o.write = TRUE)
adjacency   <- get_adjacency(mf_location)


# test that maxit works properly
test_that("Single variable GAM with network and plots", {
  # simple gam
  m_upstream  <- smnet(
    formula = Summer_mn ~ 1 + m(ELEV_DEM) + m(NEAR_X, NEAR_Y) + network(adjacency, weight = "areaPI"), 
    method = "AICC", 
    data.object = mf, 
    control = list(verbose = FALSE)
  )
  expect_equal(class(m_upstream), 'smnet')
  expect_equal(names(summary(m_upstream, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  plot(m_upstream, type = "covariates")
  plot(m_upstream, type = "segments", weight = 0.5, shadow = 3)
  plot(m_upstream, type = "full"    , weight = 0.5, shadow = 3, legend.range = c(0, 30))
  plot(m_upstream, type = "full"    , weight = 0.5, shadow = 3, legend.range = c(0, 30), sites = T)
  plot(m_upstream, res = T)
  plot(m_upstream)
  plot(m_upstream, res = T)
  plot(m_upstream, type = "segments", cex = 2)
})









