source('helpjm.r')
library(MASS)

N <- 500
b_1 <- 1
b_2 <- 2
e <- 1
e_sp <- 2
rho <- 0.01

mydat <- data.frame('xc'=runif(N, 0, 1000), 
                    'yc'=runif(N, 0, 1000),
                    'X_i'=runif(N, 0, 5))
a_1 <- rnorm(N, mean=0, sd=e)
d <- with(mydat, as.matrix(dist(cbind(xc,yc))))
a_2 <- mvrnorm(1, rep(0,N), exp(-rho*d)*e_sp^2)
mydat$y_i <- b_1 + b_2*mydat$X_i + a_1 + a_2

par(mfrow=c(1,2))
with(mydat, plot(X_i, y_i))
with(mydat, plotbubble(xc, yc, y_i, xlab='X coordinate', ylab='Y coordinate', colsign=FALSE, size=0.5))

mod1 <- lm(y_i ~ X_i, data=mydat)
round(summary(mod1)$coef,5)

library(gstat)
library(sp)
mydat$res1 <- resid(mod1)
with(mydat, plotbubble(xc, yc, res1))
coordinates(mydat) <- c('xc','yc')
with(mydat, plotvario(variogram(res1 ~ 1, data=mydat)))

sim_lm_sp <- function(N = 500,
                       b_1 = 1,
                       b_2 = 2,
                       e = 1,
                       e_sp = 2,
                       rho = 0.01) {
  mydat <- data.frame('xc'=runif(N, 0, 1000), 
                      'yc'=runif(N, 0, 1000),
                      'X_i'=runif(N, 0, 5))
  a_1 <- rnorm(N, mean=0, sd=e)
  d <- with(mydat, as.matrix(dist(cbind(xc,yc))))
  a_2 <- mvrnorm(1, rep(0,N), exp(-rho*d)*e_sp^2)
  mydat$y_i <- b_1 + b_2*mydat$X_i + a_1 + a_2
  return(mydat)
}
out <- as.data.frame(NULL)
for(i in 1:100) {
  mydat <- sim_lm_sp()
  mod <- lm(y_i ~ X_i, data=mydat)
  out <- rbind(out, c(as.vector(coef(mod)),(sum(resid(mod)^2))/(nrow(mydat)-attr(logLik(mod),'df'))))
}
par(mfrow=c(1,3))
dplot(out[,1], b_1, xlab='b_1')
dplot(out[,2], b_2, xlab='b_2')
dplot(out[,3], 1, xlab='phi')

