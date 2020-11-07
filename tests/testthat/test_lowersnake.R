context('Testing lowersnake')

mf_location <- paste0(normalizePath('testdata'), "/lsn28.ssn")
mf          <- suppressWarnings(SSN::importSSN(mf_location, predpts = NULL, o.write = TRUE))
adjacency   <- get_adjacency(mf_location)


# test that maxit works properly
test_that("Linear model with network", {
  # simple gam
  m_upstream  <- smnet(
    NEAR_Y  ~ network(adjacency, weight = "autoShreve") + NEAR_X, 
    control = list(approx = 100, verbose = FALSE), 
    data.object = mf
  )
  expect_equal(class(m_upstream), 'smnet')
  expect_equal(names(summary(m_upstream, verbose = FALSE)), c("linear.terms", "smooth.terms"))
})

