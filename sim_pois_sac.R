require(MASS)

### Simulation of POIS data WITH spatial autocorrelation:
sim_pois_sac <- function(n_obs, x_lo=0, x_hi=1, sill, range, b) {
  x <- runif(n_obs, x_lo, x_hi)
  full_grid <- sim_grid(n_obs, sill=sill, range=range)
  obs_grid <- full_grid[sample(nrow(full_grid),n_obs,replace=F),]
  dd <- data.frame(x, xc=obs_grid$x, yc=obs_grid$y)
  eta0 <- model.matrix(~x, data = dd) %*% b
  eta <- eta0 + obs_grid$sim1
  mu <- exp(eta)
  dd$y <- rpois(n_obs, lambda=mu)
  return(dd)
}

### Simulation of a series of NB data WITH spatial autocorrelation:
sim_series_pois_sac <- function(K=100, n_obs, sill, range, b) {
  ests <- as.data.frame(NULL)
  pb <- txtProgressBar(min=0, max=K, style=3)
  for (i in 1:K) {
    mydat <- sim_pois_sac(n_obs=n_obs, sill=sill, range=range, b=b)
    mod <- glm(y ~ x, data=mydat, family='poisson')
    ests <- rbind(ests, coef(mod))
    setTxtProgressBar(pb,i)
  }
  close(pb)
  return(list(ests=ests))
}
