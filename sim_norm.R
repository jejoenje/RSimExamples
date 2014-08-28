### Simulation of Gaussian data WITHOUT spatial autocorrelation:
sim_norm <- function(n_obs, x_lo=0, x_hi=1, b, s) {
  x <- runif(n_obs, x_lo, x_hi)
  # These next two lines generate x and y coordinates but note these are not used in the
  #  simulation. They are included to make output compatible with output from sim_pois_sac().
  xy <- expand.grid(xc=1:ceiling(sqrt(n_obs+1)), yc=1:ceiling(sqrt(n_obs+1)))
  xy <- xy[sample(nrow(xy),n_obs,replace=F),]
  
  dd <- data.frame(x, xy)
  eta <- model.matrix(~x, data = dd) %*% b
  mu <- eta
  dd$y <- rnorm(n_obs, mean=mu, sd=s)
  return(dd)
}
