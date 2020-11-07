context('Testing river tweed')

mf_location <- paste0(normalizePath('testdata'), "/nitrate.ssn")
mf          <- suppressWarnings(SSN::importSSN(mf_location, predpts = NULL, o.write = TRUE))
adjacency   <- get_adjacency(mf_location, netID = 2)

# extract the data frame containing pollution data inside the SSN object
data_obs          <- SSN::getSSNdata.frame(mf)
z                 <- as.Date(as.character(data_obs$Date_), "%d-%b-%y")
data_obs$year     <- data_obs$year1   <- lubridate::year(z)
data_obs$d.in.y   <- data_obs$d.in.y1 <- lubridate::yday(z)
new_tweed         <- SSN::putSSNdata.frame(data_obs, x = mf)

# test that maxit works properly
test_that("Linear model with network", {
  # simple gam
  m_upstream  <- smnet(
    formula = log(Value_) ~  m(d.in.y1, year1)+ network(adjacency, weight = "autoShreve"),
    control = list(approx = 100, 
                   maxit = 1000, 
                   optim.method = "Nelder-Mead", 
                   trace = 0, 
                   tol = 10^-8, 
                   verbose = FALSE), 
    data.object = new_tweed, 
    netID = 2)
  expect_equal(class(m_upstream), 'smnet')
  expect_equal(names(summary(m_upstream, verbose = FALSE)), c("linear.terms", "smooth.terms"))
  plot(m_upstream, type = "segments")
  plot(m_upstream, type = "covariates")
  plot(m_upstream, type = "segments", shadow = 1.5, weight = 0.8, sites = T, legend.text = "log nitrate (mg/L)")
  plot(m_upstream, type = "segments", se = T, shadow = 1.5, weight = 0.8, sites = T, legend.text = "log nitrate (mg/L)")
  
})










