
######################################################################
######################################################################
######################################################################
grid.to.ll <- function(east,north,N0=-100000,E0=400000,F0=0.9996012717,
                       phi0=pi*49/180,lambda0=-pi*2/180,
                       a=6377563.396,b=6356256.910) {
  ######################################################################
  # Function to convert grid coordinates to latitudes and longitudes.
  # Arguments as for ll.to.grid(), except that long and lat are replaced
  # by east and north respectively.
  #
  # Value: a data frame containing 2 columns named `Long' and `Lat',
  # with the translated co-ordinates of the points in east and north
  # (in degrees)
  #
  ######################################################################
  e2 <- 1 - (b/a)^2; n <- (a-b)/(a+b); M <- Inf
  phi1 <- (north-N0)/(a*F0) + phi0
  while (max(abs(north-N0-M)) >= 1e-4) {
    tmp1 <- phi1-phi0; tmp2 <- phi1+phi0
    M <- b*F0*( ((1+n+(1.25*n^2*(1+n)))*tmp1) -
                  (3*n*(1+n+(7*(n^2)/8))*sin(tmp1)*cos(tmp2)) +
                  (15*n^2*(1+n)*sin(2*tmp1)*cos(2*tmp2)/8) -
                  (35*n^3*sin(3*tmp1)*cos(3*tmp2)/24) )
    phi1 <- (north-N0-M)/(a*F0) + phi1
  }
  sp <- sin(phi1); tmp1 <- 1 - (e2*sp^2)
  nu <- a*F0/sqrt(tmp1); rho <- a*F0*(1-e2)/(tmp1^1.5)
  ratio <- (1-(e2*sp^2))/(1-e2); eta2 <- ratio-1
  T <- tan(phi1)
  z7 <- T/(2*rho*nu)
  z8 <- T*(5+(3*T^2)+eta2 -(9*eta2*T^2))/(24*rho*nu^3)
  z9 <- T*(61+(90*T^2)+(45*T^4))/(720*rho*nu^5)
  sec <- 1/cos(phi1)
  z10 <- sec/nu
  z11 <- sec*(ratio+(2*T^2))/(6*nu^3)
  z12 <- sec*(5+(28*T^2)+(24*T^4))/(120*nu^5)
  z12a <- sec*(61+(662*T^2)+(1320*T^4)+(720*T^6))/(5040*nu^7)
  tmp1 <- east-E0
  Lat <- phi1 - (z7*tmp1^2) + (z8*tmp1^4) - (z9*tmp1^6)
  Long <- lambda0 + (z10*tmp1) - (z11*tmp1^3) + (z12*tmp1^5) - (z12a*tmp1^7)
  data.frame(Long=Long*180/pi,Lat=Lat*180/pi)
}