library(MASS)

source('sim_nb.r')

K <- 100
n_obs <- 500
size <- 0.8
b <- c(1,3)

par(mfrow=c(2,2))

mydat <- sim_nb(n_obs=n_obs, size=size, b=b)
plot(y ~ x, data=mydat)

ests <- sim_series_nb(K=K, n_obs=n_obs, size=size, b=b)

apply(ests$ests, 2, mean)
mean(ests$e_size)

plot(density(ests$ests[,1]), main='b1')
lines(c(b[1],b[1]), c(0,max(density(ests$ests[,1])$y+1)), col='red', lty='dashed')
plot(density(ests$ests[,2]), main='b2')
lines(c(b[2],b[2]), c(0,max(density(ests$ests[,2])$y+1)), col='red', lty='dashed')
plot(density(ests$e_size), main='size')
lines(c(size,size), c(0,max(density(ests$e_size)$y+1)), col='red', lty='dashed')

