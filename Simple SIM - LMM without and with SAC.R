rm(list=ls())
opar <- par()
library(MASS) # for mvrnorm()
library(lme4)
library(sp)   # for bubble()
library(gstat)

source('sim_lmm_sp3.R')
source('helpjm.r')

# Parameters for simulations, definitions in sim_lmm_sp.R:
parset <-list(
  xcmin=0,
  xcmax=1000,
  ycmin=0,
  ycmax=1000,
  rho=0.01,
  re_sd=1,
  e_sp_sd=1,
  e=1,
  a=1,
  b=2,
  n_site=20,
  n_obs_site=25)

# Simulate some mixed-effects data with spatially independent response y:
parset1 <- c(parset, 'sp'=FALSE, 'spcor'=FALSE)
mydat <- do.call(sim_lmm_sp3, parset1)

# Fit model:
mod1 <- lmer(y ~ x + (1|site), data=mydat, REML=F)

#R2 test:
r2mm(mod1)
# All model estimates:
summary(mod1)
# Fixed effect estimates:
fixef(mod1)
# Random effect SD estimates:
VarCorr(mod1)
# Dispersion estimate a la Zuur:
dispZuur(mod1)
# Residuals against fitted:
par(mfrow=c(1,2))
plot(fitted(mod1),resid(mod1, type='pearson'))
# See if I can correctly calculate residuals.
# First, check predictions:
my_fit <- model.matrix(mod1) %*% fixef(mod1)
plot(my_fit, predict(mod1))
# Mmm. Nope.
my_fit <- model.matrix(mod1) %*% fixef(mod1)
# Match RE BLUP:
mydat$re1 <- ranef(mod1)$site[[1]][match(mydat$site, row.names(ranef(mod1)$site))]
my_fit <- my_fit+mydat$re1
plot(my_fit, predict(mod1)); abline(a=0, b=1, col='red')
# Bingo! So, by default, predict.lmer ADDS the BLUPs for each site, so in effect
#  the fitted predictions are specific to sites.
plot(my_fit, fitted(mod1)); abline(a=0, b=1, col='red')
# So how to work out residuals?
plot(mydat$y-fitted(mod1), resid(mod1)); abline(a=0, b=1,col='red')
# Spot on! So which is the default resid()?
plot(mydat$y-fitted(mod1), resid(mod1, type='pearson')); abline(a=0, b=1,col='red')
plot(mydat$y-fitted(mod1), resid(mod1, type='working')); abline(a=0, b=1,col='red')
plot(mydat$y-fitted(mod1), resid(mod1, type='response')); abline(a=0, b=1,col='red')
plot(mydat$y-fitted(mod1), resid(mod1, type='deviance')); abline(a=0, b=1,col='red')
# All the same??
plot(resid(mod1, type='pearson'), resid(mod1, type='working'))
plot(resid(mod1, type='pearson'), resid(mod1, type='response'))
plot(resid(mod1, type='pearson'), resid(mod1, type='deviance'))
# ...apparently...

# Do the above a bunch of times to see how well the estimation works:
#  (SAVED AS OBJECT, ONLY RUN IF NEEDED - LOAD COMMAND BELOW)

# ests <- as.data.frame(NULL)
# K <- 100
# pb <- txtProgressBar(min=0, max=K, style=3)
# for(i in 1:K) {
#   nwdat <- do.call(sim_lmm_sp3, parset1)
#   nwmod <- lmer(y ~ x + (1|site), data=nwdat, REML=F)
#   nwests <- c(as.vector(fixef(nwmod)),
#               as.data.frame(VarCorr(nwmod))[,'sdcor'],
#               dispZuur(nwmod))
#   ests <- rbind(ests, t(data.frame(nwests)))
#   setTxtProgressBar(pb,i)
# }
# close(pb)
# names(ests) <- c('a','b','re_sd','e_sd','disp')
# save(ests, file='sim_lmm_1.Rdata')

load('sim_lmm_1.Rdata')
par(mfrow=c(2,3))
dplot(ests[,'a'],parset1$a,'a')
dplot(ests[,'b'],parset1$b,'b')
dplot(ests[,'re_sd'],parset1$re_sd,'re_sd')
dplot(ests[,'e_sd'],parset1$e,'e_sd')
dplot(ests[,'disp'],1,'disp')

# Try some diag plots:
# Residuals against fitted:
plot(resid(mod1), fitted(mod1, type='response'))
# Histogram of residuals
hist(resid(mod1))
# Residuals on a map:
mydat$E <- resid(mod1)
coordinates(mydat) <- c('xc','yc')
bubble(mydat, 'E', col=c('black','grey'), main='D-residuals', xlab='X',ylab='Y')
# Variogram
mod1vg <- variogram(E ~ 1, data=mydat)
plotvario(mod1vg)

#####################
# Now consider spatially correlated data - but with RE and SAC independent:

parset2 <- c(parset, 'sp'=TRUE, 'spcor'=FALSE)
mydat2 <- do.call(sim_lmm_sp3, parset2)

# Model w/o taking SAC into account:
mod1_sp <- lmer(y ~ x + (1|site), data=mydat2, REML=F)

# All model estimates:
summary(mod1_sp)
# Fixed effect estimates:
fixef(mod1_sp)
# Random effect SD estimates:
VarCorr(mod1_sp)
# Dispersion estimate a la Zuur:
dispZuur(mod1_sp)

# Try some diag plots:
# Residuals against fitted:
plot(resid(mod1_sp), fitted(mod1_sp, type='response'))
# Histogram of residuals
hist(resid(mod1_sp))
# Residuals on a map:
mydat2$E_sp <- resid(mod1_sp)
coordinates(mydat2) <- c('xc','yc')
bubble(mydat2, 'E_sp', col=c('black','grey'), main='D-residuals', xlab='X',ylab='Y')
# Variogram
mod1_sp_vg <- variogram(E_sp ~ 1, data=mydat2)
plotvario(mod1_sp_vg)

# Run a simulation series as above:
#  (SAVED AS OBJECT, ONLY RUN IF NEEDED - LOAD COMMAND BELOW)

# ests_sp <- as.data.frame(NULL)
# K <- 100
# pb <- txtProgressBar(min=0, max=K, style=3)
# for(i in 1:K) {
#   nwdat <- do.call(sim_lmm_sp3, parset2)
#   nwmod <- lmer(y ~ x + (1|site), data=nwdat, REML=F)
#   nwests <- c(as.vector(fixef(nwmod)),
#               as.data.frame(VarCorr(nwmod))[,'sdcor'],
#               dispZuur(nwmod))
#   ests_sp <- rbind(ests_sp, t(data.frame(nwests)))
#   setTxtProgressBar(pb,i)
# }
# close(pb)
# names(ests_sp) <- c('a','b','re_sd','e_sd','disp')
# save(ests_sp, file='sim_lmm_sp_1.Rdata')

load('sim_lmm_sp_1.Rdata')
par(mfrow=c(2,3))
dplot(ests_sp[,'a'],parset2$a,'a')
dplot(ests_sp[,'b'],parset2$b,'b')
dplot(ests_sp[,'re_sd'],parset2$re_sd,'re_sd')
dplot(ests_sp[,'e_sd'],parset2$e,'e_sd')
dplot(ests_sp[,'disp'],1,'disp')




###
### Now how does clustering of sites in space affect the above.
###
rm(list=ls())
source('sim_lmm_sp.R')
source('sim_lmm_sp2.R')
source('helpjm.r')

# Parameters for simulations, definitions in sim_lmm_sp.R:
parset1 <-list(
  xcmin=0,
  xcmax=1000,
  ycmin=0,
  ycmax=1000,
  rho=0.01,
  re_sd=2,
  e_sp_sd=4,
  e=1,
  a=1,
  b=2,
  n_site=20,
  n_obs_site=100)

# Simulate some mixed-effects data with spatially independent response y and spatially correlated
#  response y_sp:
mydat2 <- do.call(sim_lmm_sp2, parset1)
plot(mydat2$xc, mydat2$yc)

# Fit LMM on SAC but not accounting for this:
mod2_sp <- lmer(y_sp ~ x + (1|site), data=mydat2, REML=F)

# All model estimates:
summary(mod2_sp)
# Fixed effect estimates:
fixef(mod2_sp)
# Random effect SD estimates:
VarCorr(mod2_sp)
# Dispersion estimate a la Zuur:
dispZuur(mod2_sp)

# Try some diag plots:
# Residuals against fitted:
plot(resid(mod2_sp), fitted(mod2_sp, type='response'))
# Histogram of residuals
hist(resid(mod2_sp))
# Residuals on a map:
mydat2$E2_sp <- resid(mod2_sp)
coordinates(mydat2) <- c('xc','yc')
bubble(mydat2, 'E2_sp', col=c('black','grey'), main='D-residuals', xlab='X',ylab='Y')
plotbubble(mydat2$xc, mydat2$yc, mydat2$E2_sp)
# example for a single site:
# site4 <- mydat2[mydat2$site==levels(mydat2$site)[4],]
# bubble(site4, 'E2_sp', col=c('black','grey'), main='D-residuals', xlab='X',ylab='Y')
# plotbubble(site4$xc, site4$yc, site4$E2_sp)
# Plot all sites:
par(mfrow=c(4,5))
par(mar=c(0.5,0.5,0.5,0.5))
for(i in 1:nlevels(mydat2$site)) {
  sitedat <- mydat2[mydat2$site==levels(mydat2$site)[i],]
  plotbubble(sitedat$xc, sitedat$yc, sitedat$E2_sp, axt='n', xlab='', ylab='')
}; par(opar)

# Variogram - I don't think this works right. Don't think these are the residuals that we need.
mod2_sp_vg <- variogram(E2_sp ~ 1, data=mydat2)
plotvario(mod2_sp_vg)


### So to check the above problem, try fitting the same model but with nlme
# First clear the decks.
rm(list=ls())
opar <- par()
library(MASS)  # for mvrnorm()
library(gstat) # for variogram()
library(sp)    # for bubble()
source('sim_lmm_sp.R')
source('sim_lmm_sp2.R')
source('helpjm.r')

# Parameters for simulations, definitions in sim_lmm_sp.R:
parset1 <-list(
  xcmin=0,
  xcmax=1000,
  ycmin=0,
  ycmax=1000,
  rho=0.01,
  re_sd=2,
  e_sp_sd=4,
  e=1,
  a=1,
  b=2,
  n_site=20,
  n_obs_site=100)

# Simulate some mixed-effects data with spatially independent response y and spatially correlated
#  response y_sp:
mydat <- do.call(sim_lmm_sp2, parset1)

# First fit the model with NLME:
library(nlme)
mod_lme <- lme(y_sp ~ x, random=~1|site, data=mydat)
mydat$E_lme <- resid(mod_lme)
plotbubble(mydat$xc, mydat$yc, mydat$E_lme)
coordinates(mydat) <- c('xc','yc')
mod_lme_vg <- variogram(E_lme ~ 1, data=mydat)
par(mfrow=c(1,1))
plotvario(mod_lme_vg)
# So the LME fit doesn't seem to give a variogram that looks spatially autocorrelated even
#  when we know the data definitely are!

# Now try the same but using lme4:
detach('package:nlme')
library(lme4)
mydat <- as.data.frame(mydat)
mod_lme4 <- lmer(y_sp ~ x + (1|site), data=mydat)
mydat$E_lme4 <- resid(mod_lme4)
plot(mydat$E_lme4, mydat$E_lme4); abline(a=0, b=1, col='red')
plotbubble(mydat$xc, mydat$yc, mydat$E_lme4)
coordinates(mydat) <- c('xc','yc')
mod_lme4_vg <- variogram(E_lme4 ~ 1, data=mydat)
par(mfrow=c(1,1))
plotvario(mod_lme4_vg)

# Clearly the variograms in the residuals are pretty similar irrespective of which model fit is used:
par(mfrow=c(1,2))
plotvario(mod_lme_vg)
plotvario(mod_lme4_vg)

###
### SO, 
# based on these simulated data, it appears that the variograms of the residuals from both LME and LME4 fits
# are very similar, irrespective of the "type" of residuals considered.
# This leaves two options - 
# 1. Both lme and lme4 residuals for some reason don't return "correct" normalised residuals as gls()
# would do.
# 2. This apparent lack of spatial autocorrelation is due to the random effect for site eating most of
# it up. In the above sims, sites are highly clustered in space.

# I could check 1. by trying to fit a GLS-type model using some extra variance structure for the RE.
#
# Alternatively, plot variograms for each site on its own, AND/OR re-do simulations where the spatial
#  and RE effects are entirely indepedent.


#### Try the latter first - re-do the sims with indepedent effects of space and RE:
### Clear the decks.
detach('package:lme4')
rm(list=ls())
opar <- par()
library(sp)   # for bubble()
source('sim_lmm_sp.R')
source('sim_lmm_sp2.R')
source('helpjm.r')
# Parameters for simulations, definitions in sim_lmm_sp.R:
parset1 <-list(
  xcmin=0,
  xcmax=1000,
  ycmin=0,
  ycmax=1000,
  rho=0.01,
  re_sd=2,
  e_sp_sd=4,
  e=1,
  a=1,
  b=2,
  n_site=20,
  n_obs_site=25)
# Simulate some mixed-effects data with spatially independent response y and spatially correlated
#  response y_sp:
mydat <- do.call(sim_lmm_sp, parset1)
# First fit the model with NLME:
library(nlme)
mod_lme <- lme(y_sp ~ x, random=~1|site, data=mydat)
mydat$E_lme <- resid(mod_lme)
plotbubble(mydat$xc, mydat$yc, mydat$E_lme)
coordinates(mydat) <- c('xc','yc')
mod_lme_vg <- variogram(E_lme ~ 1, data=mydat)
par(mfrow=c(1,1))
plotvario(mod_lme_vg)
# AHA! So a very clear pattern in the variogram in this case - suggesting clear spatial autocorrelation.
# Now repeat with lme4:
detach('package:nlme')
library(lme4)
mydat <- as.data.frame(mydat)
mod_lme4 <- lmer(y_sp ~ x + (1|site), data=mydat)
mydat$E_lme4 <- resid(mod_lme4)
plotbubble(mydat$xc, mydat$yc, mydat$E_lme4)
coordinates(mydat) <- c('xc','yc')
mod_lme4_vg <- variogram(E_lme4 ~ 1, data=mydat)
par(mfrow=c(1,1))
plotvario(mod_lme4_vg)
# Same again - as per the lme fit!
par(mfrow=c(1,2))
plotvario(mod_lme_vg)
plotvario(mod_lme4_vg)

### So it seems that a lot (or all!) of the apparent SAC in the data across all sites
###  gets completely swamped/masked by the differences between sites, at least in the Semivariogram.

### So next -  try to plot the variograms per site, instead...

### New data, fitted with lme4

rm(list=ls())
opar <- par()
library(sp)   # for bubble()
library(lme4)
source('sim_lmm_sp.R')
source('sim_lmm_sp2.R')
source('helpjm.r')
# Parameters for simulations, definitions in sim_lmm_sp.R:
parset1 <-list(
  xcmin=0,
  xcmax=1000,
  ycmin=0,
  ycmax=1000,
  rho=0.01,
  re_sd=2,
  e_sp_sd=4,
  e=1,
  a=1,
  b=2,
  n_site=20,
  n_obs_site=25)
# Simulate some mixed-effects data with spatially independent response y and spatially correlated
#  response y_sp:
mydat <- do.call(sim_lmm_sp2, parset1)
mod_lme4 <- lmer(y_sp ~ x + (1|site), data=mydat)
summary(mod_lme4)
dispZuur(mod_lme4)
mydat$E_lme4 <- resid(mod_lme4)
par(mfrow=c(1,2))
plotbubble(mydat$xc, mydat$yc, mydat$E_lme4)
coordinates(mydat) <- c('xc','yc')
mod_lme4_vg <- variogram(E_lme4 ~ 1, data=mydat)
plotvario(mod_lme4_vg)
par(mfrow=c(4,5))
par(mar=c(0.5,0.5,0.5,0.5))
for(i in 1:nlevels(mydat$site)) {
  sitedat <- mydat[mydat$site==levels(mydat$site)[i],]
  plotbubble(sitedat$xc, sitedat$yc, sitedat$E_lme4, axt='n')
}
par(mfrow=c(4,5))
par(mar=c(0.5,0.5,0.5,0.5))
for(i in 1:nlevels(mydat$site)) {
  sitedat <- mydat[mydat$site==levels(mydat$site)[i],]
  vg <- variogram(E_lme4 ~ 1, data=sitedat)
  plotvario(vg)
}

### Looks relatively strong per site - test w/o SAC:
mydat <- do.call(sim_lmm_sp2, parset1)
mod_lme4_ns <- lmer(y ~ x + (1|site), data=mydat)
summary(mod_lme4_ns)
dispZuur(mod_lme4_ns)
mydat$E_lme4_ns <- resid(mod_lme4_ns)
par(mfrow=c(1,2))
plotbubble(mydat$xc, mydat$yc, mydat$E_lme4_ns)
par(mfrow=c(4,5))
par(mar=c(0.5,0.5,0.5,0.5))
for(i in 1:nlevels(mydat$site)) {
  sitedat <- mydat[mydat$site==levels(mydat$site)[i],]
  plotbubble(sitedat$xc, sitedat$yc, sitedat$E_lme4, axt='n')
}
coordinates(mydat) <- c('xc','yc')
mod_lme4_vg_ns <- variogram(E_lme4_ns ~ 1, data=mydat)
plotvario(mod_lme4_vg_ns)
par(mfrow=c(4,5))
par(mar=c(0.5,0.5,0.5,0.5))
for(i in 1:nlevels(mydat$site)) {
  sitedat <- mydat[mydat$site==levels(mydat$site)[i],]
  vg <- variogram(E_lme4_ns ~ 1, data=sitedat)
  plotvario(vg)
}





