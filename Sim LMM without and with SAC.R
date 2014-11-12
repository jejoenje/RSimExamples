library(lme4)
library(lattice)
library(MASS)
library(ncf)

sim_lmm_sp <- function(n_obs_site=25,               # number of observations per site
                       n_sites=20,                  # number of sites (random effect levels)
                       wi_site_scale=100,           # parameter scaling the within-site variation in X and Y relative
                                                    #  to across-site variation in X and Y.
                       rho=0.01,                     # spatial corr parameter
                       re_sd=2,                     # random effect SD
                       re_sp_sd=4,                  # spatial random effect SD
                       e_sd=1,                      # residual SD
                       a=1,                         # fixed intercept
                       b=2,                         # fixed slope    
                       xmin=0,                      # minimum x predictor measurement
                       ymin=5,                       # maximum x predictor measurement
                       xcmax=382658,                # maximum x coordinate
                       xcmin=270363,                # minimum x coordinate
                       ycmax=728811,                # maximum y coordinate
                       ycmin=636634                 # minimum y coordinate
                       ) {
  
  xlims <- xcmax-xcmin
  ylims <- ycmax-ycmin
  sites <- data.frame(site=factor(paste('site',1:n_sites,sep='')), 
                      xc=round(runif(n_sites, xcmin, xcmax),0), 
                      yc=round(runif(n_sites, ycmin, ycmax),0))
  xvals <- runif(n_sites*n_obs_site, xmin, ymin)  # measured x values
  
  x <- as.vector(NULL)    
  y <- as.vector(NULL)
  site <- as.vector(NULL)
  for(i in 1:n_sites) {
    x <- c(x, round(sites[sites$site==levels(sites$site)[i],'xc']+
                      rnorm(n_obs_site, 0, xlims/wi_site_scale),0))
    y <- c(y, round(sites[sites$site==levels(sites$site)[i],'yc']+
                      rnorm(n_obs_site, 0, ylims/wi_site_scale),0))
    site <- c(site, rep(as.vector(sites$site[i]),n_obs_site))
  }
  out <- data.frame(site=site, xc=x, yc=y, xvals=xvals)
  rm(x,y,site,xvals)
  
  re <- rnorm(n_sites, 0, re_sd)
  re <- rep(re, each=n_obs_site)
  
  m <- exp(-rho*as.matrix(dist(with(out, {cbind(xc, yc)}))))        # correlation matrix for spatial RE
  re_sp <- mvrnorm(1, rep(0, n_sites*n_obs_site), m*re_sp_sd^2)     # spatial RE from multivariate N with 
                                                                    # means 0 and covariance matrix m
  
  lin_pred <- a + b*out$xvals + re + re_sp + 
              rnorm(n_sites*n_obs_site, 0, e_sd)                    # linear Y WITH spatial corr
  lin_pred_ns <- a + b*out$xvals + re + 
              rnorm(n_sites*n_obs_site, 0, e_sd)                    # linear Y WITHOUT spatial corr
  
  out$y <- lin_pred
  out$y_nospace <- lin_pred_ns
  
  return(out)  
}

  


# Fit single model to non-spatially AC data:
obsdat <- sim_lmm_sp()
plot(obsdat$xc, obsdat$yc)
mod_nospace <- lmer(y_nospace ~ xvals + (1|site), data=obsdat)
summary(mod_nospace)
# Fixed effect estimates:
fixef(mod_nospace)
# RE variance estimates:
sqrt(as.data.frame(VarCorr(mod_nospace))[,'vcov'])
# ... parameter estimates seem ok
# Plot residuals against predicted:
plot(predict(mod_nospace), resid(mod_nospace))
# Plot predicted against observed:
plot(obsdat$y_nospace, predict(mod_nospace))
# Plot 'map' of residuals:
plot(obsdat$xc, obsdat$yc, cex=resid(mod_nospace), pch=16, 
     col=gray.colors(n=length(resid(mod_nospace)))[order(resid(mod_nospace))])
# Plot 'map' of residuals for each site:
obsdat$res_nospace <- resid(mod_nospace)
par(mfrow=c(4,5))
for (i in 1:nlevels(obsdat$site)) {
  sitedat <- obsdat[obsdat$site==levels(obsdat$site)[i],]
  plot(sitedat$xc, sitedat$yc, pch=16, cex=sitedat$res_nospace, 
       col=gray.colors(n=length(sitedat$res_nospace))[order(sitedat$res_nospace)])
}
# Plot residuals against predicted for each site:
obsdat$pred_nospace <- predict(mod_nospace)
xyplot(res_nospace ~ pred_nospace | site, data=obsdat)
# Overdispersion estimate a la Zuur:
(sum(resid(mod_nospace)^2))/(nrow(obsdat)-attr(logLik(mod_nospace), "df"))

# Fit a series of models like the above:
K <- 1000
# Sim parameters for simulation series:
a <- 1
b <- 2
re_sd <- 2
e_sd <- 1
fefs <- as.data.frame(NULL)
vefs <- as.data.frame(NULL)
disp <- as.vector(NULL)
pb <- txtProgressBar(min=0, max=K, style=3)
for (i in 1:K) {
  obsdat <- sim_lmm_sp(a=a, b=b, re_sd=re_sd, e_sd=e_sd)
  mod <- lmer(y_nospace ~ xvals + (1|site), data=obsdat)
  fefs <- rbind(fefs, fixef(mod))
  vefs <- rbind(vefs, sqrt(as.data.frame(VarCorr(mod))[,'vcov']))
  disp <- c(disp, (sum(resid(mod)^2))/(nrow(obsdat)-attr(logLik(mod), "df")))
  setTxtProgressBar(pb,i)
}
close(pb)
ests <- cbind(fefs, vefs, disp)
# Estimates parameters: fixed intercept a, fixed slope b, RE variance, residual variance, dispersion:
names(ests) <- c('est_a','est_b','est_re_sd','est_e_sd','e_disp')
par(mfrow=c(2,3))
plot(density(ests$est_a), main='a'); lines(c(a,a),c(0,max(density(ests$est_a)$y)), col='red',lty='dotted')
plot(density(ests$est_b), main='b'); lines(c(b,b),c(0,max(density(ests$est_b)$y)), col='red',lty='dotted')
plot(density(ests$est_re_sd), main='re_sd'); lines(c(re_sd,re_sd),c(0,max(density(ests$est_re_sd)$y)), col='red',lty='dotted')
plot(density(ests$est_e_sd), main='e_sd'); lines(c(e_sd,e_sd),c(0,max(density(ests$est_e_sd)$y)), col='red',lty='dotted')
plot(density(ests$e_disp), main='dispersion'); lines(c(1,1),c(0,max(density(ests$e_disp)$y)), col='red',lty='dotted')
# So the parameter estimates all look pretty good.
save(ests, file='Sim LMM wo SAC.RData')



# Now do the same thing, but fit a single LMM on a single version of the data WITH autocorrelated errors,
#  but without accounting for this in the model.
# First find a value of rho (spatial correlation coefficient) that makes sense for strong correlation WI
#  sites but less so between sites.
# So what is the mean variation in X and Y coordinates within sites:
library(plyr)

xcmax<-382658
xcmin<-270363
ycmax<-728811
ycmin<-636634
wi_site_scale<-10
n_obs_site<-25
n_sites<-20
# Given the above coordinate limits and scalar vs between and within site variation...
# x <- as.vector(NULL)    
# y <- as.vector(NULL)
# site <- paste('site',1:n_sites,sep='')
# sitex <- runif(n_sites, xcmin, xcmax)
# sitey <- runif(n_sites, ycmin, ycmax)
# for(i in 1:length(site)) {
#   x <- c(x, round(rep(sitex[i],n_obs_site)+rnorm(n_obs_site, 0, (xcmax-xcmin)/wi_site_scale)))
#   y <- c(y, round(rep(sitey[i],n_obs_site)+rnorm(n_obs_site, 0, (ycmax-ycmin)/wi_site_scale)))
# }
# sites <- data.frame(site=factor(rep(as.vector(site), each=n_obs_site)), xc=x, yc=y)
# av_xydist <- mean(as.matrix(ddply(sites, .(site), summarise, 
#                                   xcdiff=max(xc)-min(xc), 
#                                   ycdiff=max(yc)-min(yc))[,2:3]))
# av_xydist
# ... av_xydist is the average distance between x and y coordinates for each site.
# So on average, there is av_xydist units of distance between observations within sites.

# The correlation between two locations is defined as a function of rho as
#  c = exp(-rho*distance)
# So if we want the correlation to be strong at close distance WITHIN sites but decrease to almost zero
#  when outside a site, 

#rho <- 4/av_xydist
#xseq <- seq(1,av_xydist, 10)
#plot(xseq, exp(-rho*xseq), xlab='Distance', ylab='Correlation')

rho <- 0.01
  
obsdat <- sim_lmm_sp(n_sites=n_sites, n_obs_site=25, rho=rho, wi_site_scale=wi_site_scale, 
                     xcmax=xcmax, xcmin=xcmin, ycmax=ycmax, ycmin=ycmin)

par(mfrow=c(2,3))
plot(obsdat$xc, obsdat$yc)
mod_space <- lmer(y ~ xvals + (1|site), data=obsdat)
summary(mod_space)
# Fixed effect estimates:
fixef(mod_space)
# RE variance estimates:
sqrt(as.data.frame(VarCorr(mod_space))[,'vcov'])
# ... parameter estimates seem ok
# Plot residuals in space:
plot(obsdat$xc, obsdat$yc, cex=resid(mod_space))
# Plot residuals against predicted:
plot(predict(mod_space), resid(mod_space))
# Plot predicted against observed:
plot(obsdat$y, predict(mod_space))
# Calculate and plot Moran's I across all data points:
mod_cl <- correlog(obsdat$xc, obsdat$xc, resid(mod_space), na.rm=T, increment=1, resamp=0)
plot(mod_cl$correlation[1:20], type='b',pch=16, cex=1.5, lwd=1.5, xlab='distance', ylab='Morans I')

# Plot residuals against predicted for each site:
obsdat$res_space <- resid(mod_space)
obsdat$pred_space <- predict(mod_space)
xyplot(res_space ~ pred_space | site, data=obsdat)
# Calculate and plot Moran's I for each site and plot:
cl_sites <- as.data.frame(NULL)
for (i in 1:nlevels(obsdat$site)) {
  sitedat <- obsdat[obsdat$site==levels(obsdat$site)[i],]
  sitedat_cl_corr <- correlog(sitedat$xc, sitedat$xc, sitedat$res_space, na.rm=T, increment=1, resamp=0)$correlation
  cl_sites <- rbind(cl_sites, data.frame(site=rep(as.vector(unique(sitedat$site)),20), 
                                         corr=sitedat_cl_corr[1:20], x=1:20))
}
xyplot(corr ~ x |site, data=cl_sites)
# Overdispersion estimate a la Zuur:
(sum(resid(mod_space)^2))/(nrow(obsdat)-attr(logLik(mod_space), "df"))



# Fit a series of models like the above:
K <- 1000
# Sim parameters for simulation series:
a <- 1
b <- 2
re_sd <- 2
e_sd <- 1
xcmax<-382658
xcmin<-270363
ycmax<-728811
ycmin<-636634
wi_site_scale<-100
n_obs_site<-25
n_sites<-20
rho<-4/av_xydist

fefs <- as.data.frame(NULL)
vefs <- as.data.frame(NULL)
disp <- as.vector(NULL)
pb <- txtProgressBar(min=0, max=K, style=3)
for (i in 1:K) {
  obsdat <- sim_lmm_sp(a=a, b=b, re_sd=re_sd, e_sd=e_sd, 
                       xcmin=xcmin, xcmax=xcmax, ycmin=ycmin, ycmax=ycmax, wi_site_scale=wi_site_scale,
                       n_obs_site, n_sites=n_sites, rho)
  mod <- lmer(y ~ xvals + (1|site), data=obsdat)
  fefs <- rbind(fefs, fixef(mod))
  vefs <- rbind(vefs, sqrt(as.data.frame(VarCorr(mod))[,'vcov']))
  disp <- c(disp, (sum(resid(mod)^2))/(nrow(obsdat)-attr(logLik(mod), "df")))
  setTxtProgressBar(pb,i)
}
close(pb)
ests_space <- cbind(fefs, vefs, disp)
# Estimates parameters: fixed intercept a, fixed slope b, RE variance, residual variance, dispersion:
names(ests_space) <- c('est_a','est_b','est_re_sd','est_e_sd','e_disp')
par(mfrow=c(2,3))
plot(density(ests_space$est_a), main='a')
lines(c(a,a),c(0,max(density(ests_space$est_a)$y)), col='red',lty='dotted')
plot(density(ests_space$est_b), main='b')
lines(c(b,b),c(0,max(density(ests_space$est_b)$y)), col='red',lty='dotted')
plot(density(ests_space$est_re_sd), main='re_sd')
lines(c(re_sd,re_sd),c(0,max(density(ests_space$est_re_sd)$y)), col='red',lty='dotted')
plot(density(ests_space$est_e_sd), main='e_sd')
lines(c(e_sd,e_sd),c(0,max(density(ests_space$est_e_sd)$y)), col='red',lty='dotted')
plot(density(ests_space$e_disp), main='dispersion')
lines(c(1,1),c(0,max(density(ests_space$e_disp)$y)), col='red',lty='dotted')
# So the parameter estimates all look pretty good.
save(ests_space, file='Sim LMM wo SAC.RData')
