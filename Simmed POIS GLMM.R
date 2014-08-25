library(lme4)
library(MASS)

sim_mm <- function(n_groups, n_obs, group_sd, b) {
  N <- n_groups*n_obs
  RE <- rnorm(n_groups, mean=0, sd=group_sd)
  x <- runif(N, 0, 1)
  dd <- data.frame(x, group=factor(paste('group',rep(1:n_groups, each=n_obs), sep='')))
  eta0 <- model.matrix(~x, data = dd) %*% b
  eta <- eta0 + RE[dd$group]
  mu <- exp(eta)
  dd$y <- rpois(N, lambda = mu)
  return(dd)
}

### Simulate K data sets and fit with LME4 Poisson:
sim_series_lme4 <- function(K=100, n_groups, n_obs, group_sd, b) {
  ests <- as.data.frame(NULL)
  pb <- txtProgressBar(min=0, max=K, style=3)
  for (i in 1:K) {
    mydat <- sim_mm(n_groups=n_groups, n_obs=n_obs, group_sd=group_sd, b=b)
    mod <- glmer(y ~ x + (1|group), data=mydat, family='poisson')
    ests <- rbind(ests, fixef(mod))
    setTxtProgressBar(pb,i)
  }
  close(pb)
  return(ests)
}

### Simulate K data sets and fit with glmmPQL:
sim_series_pql <- function(K=100, n_groups, n_obs, group_sd, b) {
  ests <- as.data.frame(NULL)
  pb <- txtProgressBar(min=0, max=K, style=3)
  for (i in 1:K) {
    mydat <- sim_mm(n_groups=n_groups, n_obs=n_obs, group_sd=group_sd, b=b)
    mod <- glmmPQL(y ~ x, random=~1|group, data=mydat, family='poisson')
    ests <- rbind(ests, fixef(mod))
    setTxtProgressBar(pb,i)
  }
  close(pb)
  return(ests)
}

b_set <- c(1,2)

### LME4 simulations and fits:
ests_lme4 <- sim_series_lme4(n_groups=10, n_obs=100, group_sd=1, b=b_set)
### Means of estimates of b
apply(ests_lme4, 2, mean)
### Mean square error of estimates of b:
apply((ests_lme4-cbind(rep(1,100),rep(2,100)))^2,2,sum)/100

### glmmPQL simulations and fits:
ests_pql <- sim_series_pql(n_groups=10, n_obs=100, group_sd=1, b=b_set)
### Means of estimates of b
apply(ests_pql, 2, mean)
### Mean square error of estimates of b:
apply((ests_pql-cbind(rep(1,100),rep(2,100)))^2,2,sum)/100




