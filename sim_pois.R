require(MASS)

### Simulation of POIS data WITHOUT spatial autocorrelation:
sim_pois <- function(n_obs, x_lo=0, x_hi=1, b) {
  x <- runif(n_obs, x_lo, x_hi)
  # These next two lines generate x and y coordinates but note these are not used in the
  #  simulation. They are included to make output compatible with output from sim_pois_sac().
  xy <- expand.grid(xc=1:ceiling(sqrt(n_obs+1)), yc=1:ceiling(sqrt(n_obs+1)))
  xy <- xy[sample(nrow(xy),n_obs,replace=F),]
  dd <- data.frame(x, xy)
  eta <- model.matrix(~x, data = dd) %*% b
  mu <- exp(eta)
  dd$y <- rpois(n_obs, lambda=mu)
  return(dd)
}

### Simulation of a series of NB data WITHOUT spatial autocorrelation:
sim_series_pois <- function(K=100, n_obs, b) {
  ests <- as.data.frame(NULL)
  pb <- txtProgressBar(min=0, max=K, style=3)
  for (i in 1:K) {
    mydat <- sim_pois(n_obs=n_obs, b=b)
    mod <- glm(y ~ x, data=mydat, family='poisson')
    ests <- rbind(ests, coef(mod))
    setTxtProgressBar(pb,i)
  }
  close(pb)
  return(list(ests=ests))
}
