source('helpjm.r')
library(MASS)
library(knitr)

###
### Common parameters:
N <- 500
b_1 <- 1
b_2 <- 2
b_3 <- 3.25
e <- 1
e_sp <- 1
rho <- 0.01
mydat <- data.frame('xc'=runif(N, 0, 1000), 
                    'yc'=runif(N, 0, 1000),
                    'X_1'=runif(N, 0, 5),
                    'X_2'=runif(N,-1,1))
a_1 <- rnorm(N, mean=0, sd=e)
d <- with(mydat, as.matrix(dist(cbind(xc,yc))))
a_2 <- mvrnorm(1, rep(0,N), exp(-rho*d)*e_sp^2)

###
### First simulate LM data with SAC:
###
mydat$y_i <- b_1 + b_2*mydat$X_1 + sin(b_3*mydat$X_2) + a_1 + a_2
par(mfrow=c(1,2))
# Plot simulated data:
with(mydat, plot(X_1, y_i))
with(mydat, plot(X_2, y_i))
with(mydat, plotbubble(xc, yc, y_i, xlab='X coordinate', ylab='Y coordinate', colsign=FALSE, size=0.5))
# Fit LM model w/o spatial structure:
mod1 <- lm(y_i ~ X_1 + X_2, data=mydat)
summary(mod1)$coef
# Plot residuals:
library(gstat)
library(sp)
mydat$res1 <- resid(mod1, type='response')
par(mfrow=c(1,2))
with(mydat, plotbubble(xc, yc, res1))
coordinates(mydat) <- c('xc','yc')
with(mydat, plotvario(variogram(res1 ~ 1, data=mydat)))
mydat <- as.data.frame(mydat)
names(mydat)[which(names(mydat)=='x')] <- 'xc'
names(mydat)[which(names(mydat)=='y')] <- 'yc'

# Fit GLS model with nlme:
library(nlme)
mod1_gls <- gls(y_i ~ X_1 + X_2, data=mydat, correlation=corExp(form=~xc+yc))  
summary(mod1_gls)  
# Plot GLS residuals:
# NOTE NEED TO USE NORMALISED RESIDUALS
mydat$res1_gls <- as.vector(resid(mod1_gls, type='normalized'))
par(mfrow=c(1,2))
plotbubble(mydat$xc, mydat$yc, mydat$res1_gls)
coordinates(mydat) <- c('xc','yc')
with(mydat, plotvario(variogram(res1_gls ~ 1, data=mydat)))
mydat <- as.data.frame(mydat)
names(mydat)[which(names(mydat)=='x')] <- 'xc'
names(mydat)[which(names(mydat)=='y')] <- 'yc'

# Fit GLM with RAC as per Crase et al 2012:
xy <- cbind(mydat$xc, mydat$yc)
library(raster)
rast <- raster(ncol=1000, nrow = 1000, ymn = 0, ymx = 1000, xmn = 0, xmx = 1000) 
res(rast) <- 1
rast[cellFromXY(rast, xy)] <- mydat$res1
focal_rac_rast <- focal(rast, matrix(1,5,5), fun=mean, na.rm=T, pad=T)
focal_rac_vect <- focal_rac_rast[cellFromXY(focal_rac_rast, xy)]
mydat <- cbind(mydat, focal_rac_vect)
mod1_rac <- lm(y_i ~ X_1 + X_2 + focal_rac_vect, data=mydat)
summary(mod1_rac)
mydat$res1_rac <- as.vector(resid(mod1_rac))
par(mfrow=c(1,2))
plotbubble(mydat$xc, mydat$yc, mydat$res1_rac)
coordinates(mydat) <- c('xc','yc')
with(mydat, plotvario(variogram(res1_rac ~ 1, data=mydat)))
mydat <- as.data.frame(mydat)
mydat <- as.data.frame(mydat)
names(mydat)[which(names(mydat)=='x')] <- 'xc'
names(mydat)[which(names(mydat)=='y')] <- 'yc'

###
### Simulate LMM data with SAC:
###
source('helpjm.r')
library(MASS)
library(knitr)
### Common parameters:
N <- 500
b_1 <- 1
b_2 <- 2
b_3 <- 3.25
e <- 1
e_sp <- 2
rho <- 0.01
mydat <- data.frame('xc'=runif(N, 0, 1000), 
                    'yc'=runif(N, 0, 1000),
                    'X_1'=runif(N, 0, 5),
                    'X_2'=runif(N, -1, 1))
a_1 <- rnorm(N, mean=0, sd=e)
d <- with(mydat, as.matrix(dist(cbind(xc,yc))))
a_2 <- mvrnorm(1, rep(0,N), exp(-rho*d)*e_sp^2)

e_g <- 1    # Between-group SD.
K <- N/25   # So the number of groups K is equal to N/25, with 25 observations per group
a_3 <- rep(rnorm(K, 0, e_g), each=25)
mydat$group <- factor(rep(paste('group',1:K,sep=''),each=25))
mydat$y2_i <- b_1 + b_2*mydat$X_1 + sin(b_3*mydat$X_2) + a_1 + a_2 + a_3
# Plot simulated data:
par(mfrow=c(2,2))
mydat <- as.data.frame(mydat)
with(mydat, plot(X_1, y2_i))
with(mydat, plot(X_2, y2_i))
with(mydat, plotbubble(xc, yc, y2_i, xlab='X coordinate', ylab='Y coordinate', colsign=FALSE, size=0.5))
# Fit LMM model:
library(lme4)
mod2 <- lmer(y2_i ~ X_1 + (1|group), data=mydat)
summary(mod2)$coefficients
# Plot residuals:
# WHICH ONES TO USE TO EVAL SAC? 'working' and 'pearson' appear to give the same.
mydat$res2 <- resid(mod2, type='pearson')
par(mfrow=c(1,2))
with(mydat, plotbubble(xc, yc, res2))
coordinates(mydat) <- c('xc','yc')
with(mydat, plotvario(variogram(res2 ~ 1, data=mydat), span=1.5))
mydat <- as.data.frame(mydat)
names(mydat)[which(names(mydat)=='x')] <- 'xc'
names(mydat)[which(names(mydat)=='y')] <- 'yc'

# Plot variograms per group:
opar <- par()
par(mfrow=c(4,5))
par(mar=c(0.1,0.1,0.1,0.1))
for(i in 1:nlevels(mydat$group)) {
  groupdat <- mydat[mydat$group==levels(mydat$group)[i],]
  coordinates(groupdat) <- c('xc','yc')
  plotvario(variogram(res2 ~ 1, data=groupdat), span=2, pts=FALSE, axes=F, xlab='',ylab='')
}
#Plot bubble plots per group:
for(i in 1:nlevels(mydat$group)) {
  groupdat <- mydat[mydat$group==levels(mydat$group)[i],]
  with(groupdat, plotbubble(xc,yc,res2,axt='n'))
}
par(opar)

# Attempt to compute RAC term, and re-fit:
xy <- cbind(mydat$xc, mydat$yc)
library(raster)
rast <- raster(ncol=1000, nrow = 1000, ymn = 0, ymx = 1000, xmn = 0, xmx = 1000) 
res(rast) <- 1
rast[cellFromXY(rast, xy)] <- mydat$res2
focal_rac_rast <- focal(rast, matrix(1,3,3), fun=mean, na.rm=T, pad=T)
focal_rac_vect <- focal_rac_rast[cellFromXY(focal_rac_rast, xy)]
mydat <- cbind(mydat, focal_rac_vect)
mod2a <- update(mod2, .~. +focal_rac_vect)
summary(mod2a)
# Plot residuals of RAC LMM:
mydat$res2a <- resid(mod2a, type='deviance')
par(mfrow=c(1,2))
with(mydat, plotbubble(xc, yc, res2a))
coordinates(mydat) <- c('xc','yc')
with(mydat, plotvario(variogram(res2a ~ 1, data=mydat)))
mydat <- as.data.frame(mydat)
names(mydat)[which(names(mydat)=='x')] <- 'xc'
names(mydat)[which(names(mydat)=='y')] <- 'yc'

### Crazy try at GLMM with extra random term instead:
mod3 <- update(mod2, .~. +(1|xc+yc))

###
### Simulate LMM data with SAC, where location correlated with group (clustered obs w/i space)
###
# Generate new, clustered coordinates for obs:
clusterscale <- 50
xc2 <- runif(nlevels(mydat$group),0,1000)
yc2 <- runif(nlevels(mydat$group),0,1000)
mydat$xc2 <- rep(xc2, each=25); rm(xc2)
mydat$yc2 <- rep(yc2, each=25); rm(yc2)
mydat$xc2 <- mydat$xc2+rnorm(nrow(mydat), 0, 1000/clusterscale)
mydat$yc2 <- mydat$yc2+rnorm(nrow(mydat), 0, 1000/clusterscale)
par(mfrow=c(1,1))
plot(mydat$xc2, mydat$yc2, xlab='X coordinate', ylab='Y coordinate')
# New distance matrix and spatial error term with new coordinates:
d2 <- with(mydat, as.matrix(dist(cbind(xc2,yc2))))
a_2 <- mvrnorm(1, rep(0,N), exp(-rho*d2)*e_sp^2)
# New response variable:
mydat$y3_i <- b_1 + b_2*mydat$X_i + a_1 + a_2 + a_3
# Plot simulated data:
par(mfrow=c(1,2))
mydat <- as.data.frame(mydat)
with(mydat, plot(X_i, y3_i))
with(mydat, plotbubble(xc2, yc2, y3_i, xlab='X coordinate', ylab='Y coordinate', colsign=FALSE, size=0.5))
# Fit LMM model:
library(lme4)
mod3 <- lmer(y3_i ~ X_i + (1|group), data=mydat)
summary(mod3)$coefficients
dispZuur(mod3)
# Plot residuals:
mydat$res3 <- resid(mod3)
par(mfrow=c(1,2))
with(mydat, plotbubble(xc2, yc2, res3))
coordinates(mydat) <- c('xc2','yc2')
with(mydat, plotvario(variogram(res3 ~ 1, data=mydat)))
mydat <- as.data.frame(mydat)
# Plot variograms per group:
opar <- par()
par(mfrow=c(4,5))
par(mar=c(0,0,0,0))
for(i in 1:nlevels(mydat$group)) {
  groupdat <- mydat[mydat$group==levels(mydat$group)[i],]
  coordinates(groupdat) <- c('xc2','yc2')
  plotvario(variogram(res3 ~ 1, data=groupdat), pts=FALSE, axes=F, xlab='',ylab='')
}

