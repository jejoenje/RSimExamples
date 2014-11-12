library(MASS)
library(ncf)    # for correlog()
library(gstat)  # for variogram()
library(sp)     # for coordinates()
library(nlme)   # for gls()

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

lin_y_sac <- a + b*mydat$x + e_sp + rnorm(N, 0, e) # Spatially correlated data

mydat$y1 <- lin_y_sac

# Fit simple GLM w/o taking SAC into account:
mod1 <- glm(y1 ~ x, data=mydat)
# Model coefficient estimates:
summary(mod1)
mydat$res1 <- resid(mod1, type='pearson')
# Manually estimate dispersion:
sum(mydat$res1^2)/(N-attr(logLik(mod1),"df"))
# Plot some diagnostics:
par(mfrow=c(2,2))
# Residuals against predicted:
plot(predict(mod1), mydat$res1)
# Plot residuals on map with shading reflecting pos/neg and size = abs size of residual:
plot(mydat$xc, mydat$yc, cex=abs(mydat$res1)/2, pch=16, 
     col=gray.colors(n=length(mydat$res1))[order(mydat$res1)])
# Calculate correlelogram and plot Moran's I:
cl_mod1 <- correlog(mydat$xc, mydat$yc, mydat$res1, increment=1, resamp=0)
plot(cl_mod1$correlation[1:20], type='o', pch=16, ylab='Moran\'s I', xlab='Distance')
#   Plot variogram OPTION 1:
# coordinates(mydat) <- c('xc','yc')
# vg1 <- variogram(res1 ~ 1, data=mydat)
# plot(vg1$dist, vg1$gamma, type='o', pch=16)
#   Plot variogram OPTION 2 - DOESN'T WORK WITH GLM?:
#vario.glm <- Variogram(mod1, form=~xc+yc, robust=T, maxDist=2000, resType='pearson')
#plot(vario.glm, smooth=T)

# Fit GLS taking SAC into account:

mod1_sac <- gls(y1 ~ x, correlation=corExp(form=~xc+yc, nugget=T), data=mydat)
summary(mod1_sac)
mydat$res1_sac <- as.vector(resid(mod1_sac, type='pearson'))
sum(mydat$res1_sac^2)/(N-attr(logLik(mod1_sac),"df"))
par(mfrow=c(2,3))
plot(predict(mod1_sac), mydat$res1_sac)
plot(mydat$xc, mydat$yc, cex=abs(mydat$res1_sac), pch=16, 
     col=gray.colors(n=length(mydat$res1_sac))[order(mydat$res1_sac)])
cl_mod1_sac <- correlog(mydat$xc, mydat$yc, mydat$res1_sac, increment=1, resamp=0)
plot(cl_mod1_sac$correlation[1:20], type='o', pch=16, ylab='Moran\'s I', xlab='Distance')
vg1_sac <- variogram(res1_sac ~ 1, data=mydat)
plot(vg1_sac$dist, vg1_sac$gamma, type='o', pch=16)
vario.GLS <- Variogram(mod1_sac, form=~xc+yc, robust=T, maxDist=2000, resType='pearson')
plot(vario.GLS, smooth=T)
