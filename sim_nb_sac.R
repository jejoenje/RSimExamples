require(MASS)

### Simulation of a NB WITH spatial autocorrelation.
### Sill and range parameters as in sim_grid() function.
sim_nb_sac <- function(n_obs, sill=sill, range=range, size, b) {
  x <- runif(n_obs, 0, 1)
  full_grid <- sim_grid(n_obs, sill=sill, range=range)
  obs_grid <- full_grid[sample(nrow(full_grid),n_obs,replace=F),]
  dd <- data.frame(x=x, xc=obs_grid$x, yc=obs_grid$y)
  eta0 <- model.matrix(~x, data = dd) %*% b
  eta <- eta0 + obs_grid$sim1
  mu <- exp(eta)
  dd$y <- rnbinom(n_obs, size = size, prob=size/(size+mu))
  return(dd)
}

### Simulation of a series of NB WITH spatial autocorrelation:
### Sill and range parameters as in sim_grid() function.
sim_series_nb_sac <- function(K=100, sill=sill, range=range, n_obs, size, b) {
  ests <- as.data.frame(NULL)
  e_size <- as.vector(NULL)
  pb <- txtProgressBar(min=0, max=K, style=3)
  for (i in 1:K) {
    mydat <- sim_nb_sac(n_obs=n_obs, sill=sill, range=range, size=size, b=b)
    mod <- glm.nb(y ~ x, data=mydat)
    ests <- rbind(ests, coef(mod))
    e_size <- c(e_size, mod$theta)
    setTxtProgressBar(pb,i)
  }
  close(pb)
  return(list(ests=ests, e_size=e_size))
}