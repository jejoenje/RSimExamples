### Simulate data Y from Poisson distribution as function of x,
###  fit GLM

source('sim_grid.r')
source('sim_pois.r')

K <- 100
n_obs <- 500
b <- c(-1,1)
x_lo <- 2
x_hi <- 4

par(mfrow=c(1,1))
mydat <- sim_pois(n_obs=n_obs, x_lo=x_lo, x_hi=x_hi, b=b)
plot(y ~ x, data=mydat)
mod <- glm(y ~ x, data=mydat, family='poisson')
nwx <- data.frame(x=seq(x_lo,x_hi,(x_hi-x_lo)/100))
pred <- predict(mod, newdata=nwx, type='response', se=T)
lines(nwx$x, pred$fit, lwd=2, col='red')
lines(nwx$x, pred$fit+2*pred$se, lwd=1, col='red')
lines(nwx$x, pred$fit-2*pred$se, lwd=1, col='red')
ci <- confint(mod)
lines(nwx$x, exp(model.matrix(~x, data=nwx) %*% ci[,2]), lty='dashed', col='red')
lines(nwx$x, exp(model.matrix(~x, data=nwx) %*% ci[,1]), lty='dashed', col='red')

par(mfrow=c(2,3))

mydat <- sim_pois(n_obs=n_obs, x_lo=x_lo, x_hi=x_hi, b=b)
plot(y ~ x, data=mydat)

mod <- glm(y ~ x, data=mydat, family='poisson')
hist(resid(mod))
plot(predict(mod, type='link'), resid(mod))

ests <- sim_series_pois(K=K, n_obs=n_obs, b=b)

apply(ests$ests, 2, mean)

plot(density(ests$ests[,1]), main='b1')
lines(c(b[1],b[1]), c(0,max(density(ests$ests[,1])$y+1)), col='red', lty='dashed')
plot(density(ests$ests[,2]), main='b2')
lines(c(b[2],b[2]), c(0,max(density(ests$ests[,2])$y+1)), col='red', lty='dashed')

