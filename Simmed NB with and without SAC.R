library(MASS)
library(ncf)
source('sim_grid.r')
source('sim_nb.r')
source('sim_nb_sac.r')

K <- 100
n_obs <- 1000
sill <- 1
range <- 10
size <- 0.8
b <- c(1,3)


### Single simulation of NONE-SAC NB:
par(mfrow=c(2,3))
mydat <- sim_nb(n_obs=n_obs, size=size, b=b)
plot(mydat$x, mydat$y)
mod <- glm.nb(y ~ x, data=mydat); coef(mod)
hist(resid(mod))
plot(predict(mod, type='link'), resid(mod))
plot(mydat$xc, mydat$yc, cex=resid(mod))
mod_corr <- correlog(mydat$xc, mydat$yc, resid(mod), na.rm=T, increment=1, resamp=0)
plot(mod_corr$correlation[1:20], type='o', pch=16, cex=1, lwd=1.5, 
     xlab='distance', ylab='Morans I', cex.lab=1,cex.axis=1)
coordinates(mydat) <- c('xc','yc')
mydat_vario <- variogram(resid(mod)~1, data=mydat)
plot(gamma ~ dist, data=mydat_vario, ylim=c(0,1.1), pch=16, type='o')

### Single simulation of SAC NB:
par(mfrow=c(2,3))
mydat_sac <- sim_nb_sac(n_obs=n_obs, sill=sill, range=range, size=size, b=b)
mod_sac <- glm.nb(y ~ x, data=mydat_sac); coef(mod_sac)
plot(mydat_sac$x, mydat_sac$y)
hist(resid(mod_sac))
plot(predict(mod_sac, type='link'), resid(mod_sac))
plot(mydat_sac$xc, mydat_sac$yc, cex=resid(mod_sac))
mod_sac_corr <- correlog(mydat_sac$xc, mydat_sac$yc, resid(mod_sac), 
                         na.rm=T, increment=1, resamp=0)
plot(mod_sac_corr$correlation[1:20], type='o', pch=16, cex=1, lwd=1.5, 
     xlab='distance', ylab='Morans I', cex.lab=1,cex.axis=1)
coordinates(mydat_sac) <- c('xc','yc')
mydat_sac_vario <- variogram(resid(mod_sac)~1, data=mydat_sac)
plot(gamma ~ dist, data=mydat_sac_vario, ylim=c(0,1.1), pch=16, type='o')


### Series of NON-SAC NB.
### Plot estimated distributions with actual parameters.
par(mfrow=c(2,3))
nb_series <- sim_series_nb(K=K, n_obs=n_obs, size=size, b=b)
plot(density(nb_series$ests[,1]))
lines(c(b[1],b[1]), c(0,max(density(nb_series$ests[,1])$y+1)), col='red', lty='dashed')
plot(density(nb_series$ests[,2]))
lines(c(b[2],b[2]), c(0,max(density(nb_series$ests[,2])$y+1)), col='red', lty='dashed')
plot(density(nb_series$e_size))
lines(c(size,size), c(0,max(density(nb_series$e_size)$y+1)), col='red', lty='dashed')

### Series of SAC NB.
### Plot estimated distributions with actual parameters.
nbsac_series <- sim_series_nb_sac(K=K, n_obs=n_obs, sill=sill, range=range, size=size, b=b)
plot(density(nbsac_series$ests[,1]))
lines(c(b[1],b[1]), c(0,max(density(nbsac_series$ests[,1])$y+1)), col='red', lty='dashed')
plot(density(nbsac_series$ests[,2]))
lines(c(b[2],b[2]), c(0,max(density(nbsac_series$ests[,2])$y+1)), col='red', lty='dashed')
plot(density(nbsac_series$e_size))
lines(c(size,size), c(0,max(density(nbsac_series$e_size)$y+1)), col='red', lty='dashed')
