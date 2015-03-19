require(MASS)

### Simulation of NB data WITHOUT spatial autocorrelation:
sim_nb <- function(n_obs, size, b) {
  x <- runif(n_obs, 0, 1)
  # These next two lines generate x and y coordinates but note these are not used in the
  #  simulation. They are included to make output compatible with output from sim_nb_sac().
  xy <- expand.grid(xc=1:ceiling(sqrt(n_obs+1)), yc=1:ceiling(sqrt(n_obs+1)))
  xy <- xy[sample(nrow(xy),n_obs,replace=F),]
  dd <- data.frame(x, xy)
  eta <- model.matrix(~x, data = dd) %*% b
  mu <- exp(eta)
  dd$y <- rnbinom(n_obs, size = size, prob=size/(size+mu))
  return(dd)
}

### Simulation of a series of NB data WITHOUT spatial autocorrelation:
sim_series_nb <- function(K=100, n_obs, size, b) {
  ests <- as.data.frame(NULL)
  e_size <- as.vector(NULL)
  pb <- txtProgressBar(min=0, max=K, style=3)
  for (i in 1:K) {
    mydat <- sim_nb(n_obs=n_obs, size=size, b=b)
    mod <- glm.nb(y ~ x, data=mydat)
    ests <- rbind(ests, coef(mod))
    e_size <- c(e_size, mod$theta)
    setTxtProgressBar(pb,i)
  }
  close(pb)
  return(list(ests=ests, e_size=e_size))
}
