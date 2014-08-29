source('sim_norm.r')

n_obs <- 500
x_lo <- 1
x_hi <- 2
b <- c(15,2)
s <- 3

mydat <- sim_norm(n_obs, x_lo=x_lo, x_hi=x_hi, b=b, s=s)
plot(mydat$x, mydat$y)

summary(mod <- glm(y ~ x, data=mydat))

nwx <- data.frame(x=seq(x_lo, x_hi, (x_hi-x_lo)/100))

### Predictions from predict():
pred1 <- predict(mod, newdata=nwx, se=T)
lines(nwx$x, pred1$fit, lwd=2, col='red')
lines(nwx$x, pred1$fit-pred1$se*1.96, col='red', lty='dashed')
lines(nwx$x, pred1$fit+pred1$se*1.96, col='red', lty='dashed')
pred2 <- predict(mod, newdata=nwx, se=T, type='response')
lines(nwx$x, pred2$fit, lwd=2, col='blue')
lines(nwx$x, pred2$fit-pred1$se*1.96, col='blue', lty='dashed')
lines(nwx$x, pred2$fit+pred1$se*1.96, col='blue', lty='dashed')

### Predictions from confint() upper and lower levels:
coef_ci <- confint(mod)
nwx_mm <- model.matrix(~x, data=nwx)
lines(nwx$x, nwx_mm %*% coef_ci[,1], col='darkgreen', lty='dashed')
lines(nwx$x, nwx_mm %*% coef_ci[,2], col='darkgreen', lty='dashed')

### Predictions from ARM sim():
library(arm)
mod_sim <- sim(mod, 1000)
coef_sim <- coef(mod_sim)
pred_sim <- as.data.frame(NULL)
for (i in 1:nrow(coef_sim)) {
  pred_sim <- rbind(pred_sim, t(nwx_mm %*% coef_sim[i,]))
}
pred_sim_mn <- apply(pred_sim, 2, mean)
pred_sim_lo <- apply(pred_sim, 2, function(x) quantile(x, probs=c(0.025)))
pred_sim_hi <- apply(pred_sim, 2, function(x) quantile(x, probs=c(0.975)))
lines(nwx$x, pred_sim_mn, lwd=2, col='purple')
lines(nwx$x, pred_sim_lo, lwd=1, col='purple', lty='dotted')
lines(nwx$x, pred_sim_hi, lwd=1, col='purple', lty='dotted')

### Predictions from home-brew sim() equivalent (to check) with sim() output):
library(mvtnorm)

pred_simV1 <- function(fit, nsims=1000, xv) {
  
  fit_df <- attr(logLik(fit),'df')      # = n-k = estimated number of parameters.
  d <- summary(fit)$dispersion          # sigma_hat^2; estimated residual VARIANCE
  xmn <- min(model.matrix(fit)[,xv])
  xmx <- max(model.matrix(fit)[,xv])
  xdf <- data.frame(x=seq(xmn, xmx, (xmx-xmn)/100))
  xmm <- model.matrix(~x, data=xdf)
  eta <- as.data.frame(NULL)
  for(i in 1:nrow(xmm)) {
    x_j <- as.vector(NULL)
    for(j in 1:nsims) {
      coefs <- rmvnorm(1, coef(fit), vcov(fit))
      x_j <- c(x_j, (xmm[i,] %*% t(coefs)))
    }
    eta <- rbind(eta, x_j)
  }
  return(t(eta))
}

pred_simV2 <- function(fit, nsims=1000, xv) {
  
  fit_df <- attr(logLik(fit),'df')      # = n-k = estimated number of parameters.
  d <- summary(fit)$dispersion          # sigma_hat^2; estimated residual VARIANCE
  s <- sqrt(d)*sqrt((fit_df)/rchisq(nsims, fit_df)) # sigma; simulated residual SD
  xmn <- min(model.matrix(fit)[,xv])
  xmx <- max(model.matrix(fit)[,xv])
  xdf <- data.frame(x=seq(xmn, xmx, (xmx-xmn)/100))
  xmm <- model.matrix(~x, data=xdf)
  eta <- as.data.frame(NULL)
  for(i in 1:nrow(xmm)) {
    x_j <- as.vector(NULL)
    for(j in 1:nsims) {
      coefs <- rmvnorm(1, coef(fit), vcov(fit)*s[j]^2)
      x_j <- c(x_j, (xmm[i,] %*% t(coefs)))
    }
    eta <- rbind(eta, x_j)
  }
  return(t(eta))
}

pred_jm1 <- pred_simV1(mod, nsims=1000, 'x')
pred_jm1_mn <- apply(pred_jm1, 2, mean)
pred_jm1_lo <- apply(pred_jm1, 2, function(x) quantile(x, probs=c(0.025)))
pred_jm1_hi <- apply(pred_jm1, 2, function(x) quantile(x, probs=c(0.975)))
lines(nwx$x, pred_jm1_mn, lwd=2, col='blue')
lines(nwx$x, pred_jm1_lo, lwd=1, col='blue', lty='dotted')
lines(nwx$x, pred_jm1_hi, lwd=1, col='blue', lty='dotted')

pred_jm2 <- pred_simV2(mod, nsims=1000, 'x')
pred_jm2_mn <- apply(pred_jm2, 2, mean)
pred_jm2_lo <- apply(pred_jm2, 2, function(x) quantile(x, probs=c(0.025)))
pred_jm2_hi <- apply(pred_jm2, 2, function(x) quantile(x, probs=c(0.975)))
lines(nwx$x, pred_jm2_mn, lwd=2, col='darkgreen')
lines(nwx$x, pred_jm2_lo, lwd=1, col='darkgreen', lty='dotted')
lines(nwx$x, pred_jm2_hi, lwd=1, col='darkgreen', lty='dotted')
## overplot upper and lower levels as produced by sim():
lines(nwx$x, pred_sim_lo, lwd=1, col='purple', lty='dotted')
lines(nwx$x, pred_sim_hi, lwd=1, col='purple', lty='dotted')

## Directly compare sim() uppper and lower levels and the predictions from pred_simV1():
plot(pred_sim_lo, pred_jm1_lo)
points(pred_sim_hi, pred_jm1_hi, col='red')
plot(nwx$x, pred_sim_lo, type='l', ylim=c(min(c(pred_sim_lo,pred_jm1_lo)), max(c(pred_sim_hi,pred_jm1_hi))))
lines(nwx$x, pred_sim_hi, type='l')
lines(nwx$x, pred_jm1_lo, type='l', col='red')
lines(nwx$x, pred_jm1_hi, type='l', col='red')

### So, pred_simV1() produces -nearly- identical predictions to using sim().
### Interestingly this version of the algorithm IGNORES variability in the estimate of
###  sigma - so the covariance matrix used to simulate coefficient values from the MVR
###  is not scaled by the extent of the "error variance".
### Regardless, predicted lines are "wobbly" instead of the smooth result produced by
###  predictions from sim(). Why?

### pred_simV1a() is the same as pred_simV1() but returns coefficient estimates rather
###  than predictions.
coef_simV1a <- function(fit, nsims=1000, xv) {
  ### Tested 28/08/2014 for Gaussian models - works as sim()
  ### Tested 28/08/2014 for Poisson models - works as sim()
  coefs <- as.data.frame(NULL)
  for(j in 1:nsims) {
    coefs <- rbind(coefs, rmvnorm(1, coef(fit), vcov(fit)))
  }
  return(coefs)
}

coef_simV1b <- function(fit, nsims=1000, xv) {
  fit_df <- attr(logLik(fit),'df')      # = n-k = estimated number of parameters.
  d <- summary(fit)$dispersion          # sigma_hat^2; estimated residual VARIANCE
  s <- sqrt(d)*sqrt((fit_df)/rchisq(nsims, fit_df)) # sigma; simulated residual SD
  coefs <- as.data.frame(NULL)
  for(j in 1:nsims) {
    coefs <- rbind(coefs, rmvnorm(1, coef(fit), vcov(fit)*s[j]^2))
  }
  return(coefs)
}

### Simulate coefficients using the homebrew algorithms above _jm3 as before (ignoring estimated
###  residual variance) and _jm4 including it.
### Also simulate using sim(), and compare distributions of estimates using density plots:
coef_jm3 <- coef_simV1a(mod, nsims=1000, 'x')
coef_jm4 <- coef_simV1b(mod, nsims=1000, 'x')
coef_sim <- coef(sim(mod, 1000))
dens_jm3_b1 <- density(coef_jm3[,1])
dens_jm4_b1 <- density(coef_jm4[,1])
dens_sim_b1 <- density(coef_sim[,1])
dens_jm3_b2 <- density(coef_jm3[,2])
dens_jm4_b2 <- density(coef_jm4[,2])
dens_sim_b2 <- density(coef_sim[,2])
dens_b1_mx <- max(c(dens_jm3_b1$y, dens_jm4_b1$y, dens_sim_b1$y))
dens_b2_mx <- max(c(dens_jm3_b2$y, dens_jm4_b2$y, dens_sim_b2$y))
par(mfrow=c(1,2))
plot(dens_sim_b1, ylim=c(0,dens_b1_mx))
lines(dens_jm3_b1, col='red')
lines(dens_jm4_b1, col='purple')
plot(dens_sim_b2, ylim=c(0,dens_b2_mx))
lines(dens_jm3_b2, col='red')
lines(dens_jm4_b2, col='purple')
### So it is clear that coef_simV1a() does the same thing as sim(), which means sim() ignores the estimated 
###  error variance.

### So can we get the same predictions from sim() and coef_simV1a()?

pred_sim <- as.data.frame(NULL)
pred_jm3 <- as.data.frame(NULL)
for (i in 1:nrow(coef_sim)) {
  pred_sim <- rbind(pred_sim, t(nwx_mm %*% coef_sim[i,]))
  pred_jm3 <- rbind(pred_jm3, t(nwx_mm %*% as.numeric(coef_jm3[i,])))
  print(i)
}
par(mfrow=c(1,1))
plot(mydat$x, mydat$y)
lines(nwx$x, apply(pred_sim, 2, mean), lwd=2, col='black')
lines(nwx$x, apply(pred_sim, 2, function(x) quantile(x, probs=c(0.025))), lwd=1.5, col='black', lty='dashed')
lines(nwx$x, apply(pred_sim, 2, function(x) quantile(x, probs=c(0.975))), lwd=1.5, col='black', lty='dashed')
lines(nwx$x, apply(pred_jm3, 2, mean), lwd=2, col='red')
lines(nwx$x, apply(pred_jm3, 2, function(x) quantile(x, probs=c(0.025))), lwd=1.5, col='red', lty='dashed')
lines(nwx$x, apply(pred_jm3, 2, function(x) quantile(x, probs=c(0.975))), lwd=1.5, col='red', lty='dashed')

### Bingo - so coef_simV1a() exactly replicates the behaviour of sim().

### So now try to do the same thing with the coefficients simulated taking into account
###  residual variance (pred_simV1b())
coef_jm4 <- coef_simV1b(mod, nsims=5000, 'x')
pred_jm4 <- as.data.frame(NULL)
for (i in 1:nrow(coef_jm4)) {
  pred_jm4 <- rbind(pred_jm4, t(nwx_mm %*% as.numeric(coef_jm4[i,])))
  print(i)
}
plot(mydat$x, mydat$y)
lines(nwx$x, apply(pred_jm4, 2, mean), lwd=2, col='red')
lines(nwx$x, apply(pred_jm4, 2, function(x) quantile(x, probs=c(0.025))), lwd=1, col='red')
lines(nwx$x, apply(pred_jm4, 2, function(x) quantile(x, probs=c(0.975))), lwd=1, col='red')

### More "wriggly" when simulating fewer cases:
plot(mydat$x, mydat$y)
coef_jm5 <- coef_simV1b(mod, nsims=1000, 'x')
pred_jm5 <- as.data.frame(NULL)
for (i in 1:nrow(coef_jm5)) {
  pred_jm5 <- rbind(pred_jm5, t(nwx_mm %*% as.numeric(coef_jm5[i,])))
  print(i)
}
lines(nwx$x, apply(pred_jm5, 2, mean), lwd=1.5, col='orange')
lines(nwx$x, apply(pred_jm5, 2, function(x) quantile(x, probs=c(0.025))), lwd=1, col='orange', lty='dotted')
lines(nwx$x, apply(pred_jm5, 2, function(x) quantile(x, probs=c(0.975))), lwd=1, col='orange', lty='dotted')

lines(loess.smooth(nwx$x, apply(pred_jm5, 2, function(x) quantile(x, probs=c(0.025))))$x,
  loess.smooth(nwx$x, apply(pred_jm5, 2, function(x) quantile(x, probs=c(0.025))))$y, col='blue')
lines(loess.smooth(nwx$x, apply(pred_jm5, 2, function(x) quantile(x, probs=c(0.975))))$x,
      loess.smooth(nwx$x, apply(pred_jm5, 2, function(x) quantile(x, probs=c(0.9755))))$y, col='blue')


### FINALLY, AN ALTERNATIVE APPROACH
### 
### FROM http://glmm.wikidot.com/faq
# Predictions and/or confidence (or prediction) intervals on predictions
# 
# The general recipe for computing predictions from a linear or generalized linear model is to:
# 
# -figure out the model matrix X corresponding to the new data;
# -* matrix-multiply X by the parameter vector β to get the predictions (or linear predictor in the case of GLM(M)s);
# -extract the variance-covariance matrix of the parameters V
# -compute XVX′ to get the variance-covariance matrix of the predictions;
# -extract the diagonal of this matrix to get variances of predictions;
# -if computing prediction rather than confidence intervals, add the residual variance;
# -take the square-root of the variances to get the standard deviations (errors) of the predictions;
# -compute confidence intervals based on a Normal approximation;
# -for GL(M)Ms, run the confidence interval boundaries (not the standard errors) through the inverse-link function.
# 

eta0 <- nwx_mm %*% coef(mod)
mod_vcov <- vcov(mod)
eta0_vcov <- nwx_mm %*% mod_vcov %*% t(nwx_mm)
eta0_var <- diag(eta0_vcov)+summary(mod)$dispersion
eta0_sd <- sqrt(eta0_var)
eta_lo <- eta0-eta0_sd
eta_hi <- eta0+eta0_sd
eta_mu <- eta0

lines(nwx$x, eta_mu, lwd=2, col='blue', lty='dashed')
lines(nwx$x, eta_hi, lwd=1, col='blue', lty='dotted')
lines(nwx$x, eta_lo, lwd=1, col='blue', lty='dotted')

### Add confint() based profile likelihood CI's for comparison:
coef_ci <- confint(mod)
lines(nwx$x, nwx_mm %*% coef_ci[,1], col='darkgreen', lty='dashed')
lines(nwx$x, nwx_mm %*% coef_ci[,2], col='darkgreen', lty='dashed')

