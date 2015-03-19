### Simulate data Y from Poisson distribution as function of x,
###  fit GLM

source('sim_grid.r')
source('sim_pois_sac.r')

K <- 100
n_obs <- 500
sill <- 0.025
range <- 5
b <- c(-1,2)
x_lo <- 1
x_hi <- 6

### Sigle simulation and model fit to show residual patterns etc.
par(mfrow=c(3,3))
# Generate data
mydat <- sim_pois_sac(n_obs=n_obs, x_lo=1, x_hi=2, sill=sill, range=range, b=b)
# Plot simulated measured data (response) y against sole explanatory variable x:
plot(y ~ x, data=mydat)
# Fit GLM with Poisson distribution to simulated data:
mod <- glm(y ~ x, data=mydat, family='poisson')
# Overplot fitted predictions:
nwx <- data.frame(x=seq(x_lo,x_hi,0.01))
pred <- predict(mod, type='response', newdata=nwx, se=T)
lines(nwx$x, pred$fit, lwd=2, col='red')
lines(nwx$x, pred$fit+2*pred$se, col='red')
lines(nwx$x, pred$fit-2*pred$se, col='red')
# Histogram of the residuals:
hist(resid(mod))
# Residuals against fitted (check for homogeneity)
plot(predict(mod, type='link'), resid(mod))
# Plot residuals on x/y coordinate system. This should show clear clustering.
plot(mydat$xc, mydat$yc, cex=resid(mod))
# Calculate correlelogram and plot Moran's I, as a visual check for spatial
#  autocorrelation:
mod_corr <- correlog(mydat$xc, mydat$yc, resid(mod), 
                         na.rm=T, increment=1, resamp=0)
plot(mod_corr$correlation[1:20], type='o', pch=16, cex=1, lwd=1.5, 
     xlab='distance', ylab='Morans I', cex.lab=1,cex.axis=1)
# In addition to Moran's I, make a semi-variogram: 
coordinates(mydat) <- c('xc','yc')
mydat_vario <- variogram(resid(mod)~1, data=mydat)
plot(gamma ~ dist, data=mydat_vario, pch=16, type='o')

### Now repeat above data simulation and model fit K times, each time
###  extract estimated parameter values.
ests <- sim_series_pois_sac(K=K, n_obs=n_obs, sill=sill, range=range, b=b)
# Calculate mean estimates over K rungs:
apply(ests$ests, 2, mean)
# Show the distribution of estimated parameter values as densityplots and add
#  a line for the position of the 'true' estimate:
plot(density(ests$ests[,1]), main='b1')
lines(c(b[1],b[1]), c(0,max(density(ests$ests[,1])$y+1)), col='red', lty='dashed')
plot(density(ests$ests[,2]), main='b2')
lines(c(b[2],b[2]), c(0,max(density(ests$ests[,2])$y+1)), col='red', lty='dashed')



