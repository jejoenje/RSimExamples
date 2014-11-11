library(glmmADMB)
library(lattice)

sim_mm_nb <- function(n_groups, n_obs, group_sd, size, b) {
  N <- n_groups*n_obs
  RE <- rnorm(n_groups, mean=0, sd=group_sd)
  x <- runif(N, 0, 1)
  dd <- data.frame(x, group=factor(paste('group',rep(1:n_groups, each=n_obs), sep='')))
  eta0 <- model.matrix(~x, data = dd) %*% b
  eta <- eta0 + RE[dd$group]
  mu <- exp(eta)
  dd$y <- rnbinom(N, size = size, prob=size/(size+mu))
  return(dd)
}

sim_series_admb <- function(K=100, n_groups, n_obs, group_sd, size, b) {
  ests <- as.data.frame(NULL)
  e_var <- as.vector(NULL)
  e_size <- as.vector(NULL)
  pb <- txtProgressBar(min=0, max=K, style=3)
  for (i in 1:K) {
    mydat <- sim_mm_nb(n_groups=n_groups, n_obs=n_obs, group_sd=group_sd, size=size, b=b)
    mod <- glmmadmb(y ~ x + (1|group), data=mydat, family='nbinom')
    ests <- rbind(ests, fixef(mod))
    e_var <- c(e_var, as.numeric(sqrt(VarCorr(mod)[[1]])))
    e_size <- c(e_size, mod$alpha)
    setTxtProgressBar(pb,i)
  }
  close(pb)
  return(list(ests=ests, e_var=e_var, e_size=e_size))
}

# single model
set.seed(52358)
mydat <- sim_mm_nb(n_groups=12,n_obs=100, group_sd=1, size=0.8, b=c(1,2))
xyplot(y ~ x | group, data=mydat)
mod <- glmmadmb(y ~ x + (1|group), data=mydat, family='nbinom')

# series of sims and fits
ests <- sim_series_admb(n_groups=12, n_obs=50, group_sd=1, size=0.8, b=c(1,2))

save(ests, 'sim_series_mmnb_admb.Rdata')
apply(ests$ests, 2, mean)
mean(ests$e_var)
mean(ests$e_size)


