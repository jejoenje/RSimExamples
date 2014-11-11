library(ncf)    # for correlog()
library(gstat)  # for variogram()
library(sp)     # for coordinates()

N <- 500
xcmin <- 0
xcmax <- 1000
ycmin <- 0
ycmax <- 1000
rho <- 0.01
e_sp_sd <- 4
e <- 1
a <- 1
b <- 2

mydat <- data.frame(xc=runif(N, xcmin, xcmax), yc=runif(N, ycmin, ycmax), x=runif(N, 0,5))

d <- as.matrix(dist(cbind(mydat$xc, mydat$yc)))
m <- exp(-rho*d)
e_sp <- mvrnorm(1, rep(0, N), m*e_sp_sd^2)

lin_y <- a + b*mydat$x + e_sp + rnorm(N, 0, e) # Spatially correlated data
#lin_y <- a + b*mydat$x + rnorm(N, 0, e)        # Spatially independent data

mydat$y1 <- lin_y

mod1 <- glm(y1 ~ x, data=mydat)
summary(mod1)
sum(resid(mod1)^2)/(N-2)
par(mfrow=c(2,2))
plot(predict(mod1), rstandard(mod1))
plot(mydat$xc, mydat$yc, cex=rstandard(mod1), pch=16, col=gray.colors(n=length(rstandard(mod1)))[order(rstandard(mod1))])
plot(mydat$y1, predict(mod1))
cl_mod1 <- correlog(mydat$xc, mydat$yc, rstandard(mod1), increment=1, resamp=0)
plot(cl_mod1$correlation[1:20], type='b')


coordinates(mydat) <- c('xc','yc')
plot(variogram(rstandard(mod1) ~ 1, data=mydat))
