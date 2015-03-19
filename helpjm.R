#helpjm.r

# JM HELPER FUNCTIONS

is.factor.df <- function(x) {
  out <- as.vector(NULL)
  for(i in 1:ncol(x)) {
    out <- c(out, is.factor(x[,i]))
  }
  return(out)
}

# Function to plot density curve of given vector dvect,
#  adding a vertical line at some reference point pval and naming the
#  plot 'pname'.
dplot <- function(dvect, pval, pname='', xlab=NULL) {
  d <- density(dvect)
  if(is.null(xlab)) {
    plot(d, main=pname)
  } else {
    plot(d, main=pname, xlab=xlab)
  }
  dmax <- max(d$y)
  lines(c(pval,pval),c(0,dmax),col='red',lty='dotted')
}

# Calculate overdispersion estimate a la Zuur 2009 (pp 224), from glm, glmer fit.
# This uses the ratio between the sum of the squared Pearson residuals (which is 
#   equivalent to the Pearson Chi-squared goodness of fit statistic) and the residual 
#   degrees of freedom, and should ideally be near 1.
# Note that this is a rough estimate only as arguably this underestimates the number
#   of "parameters" inferred by the random structure of any model.
dispZuur <- function(mod, restype='pearson') {
  if(restype=='pearson' | restype=='deviance') {
    if(restype=='pearson') {
      return(sum(resid(mod, type='pearson')^2)/(nrow(model.matrix(mod))-attr(logLik(mod),'df')))  
    }
    if(restype=='deviance') {
      return(sum(resid(mod, type='deviance')^2)/(nrow(model.matrix(mod))-attr(logLik(mod),'df')))  
    }
  } else {
    print('Error - residual type should be either \'pearson\' or \'deviance\'')
  }
}

r2mm <- function(mod) {
  if(is(mod,'glmerMod')) {
    mod_fam <- family(mod)$family
    mod_lnk <- family(mod)$link
  }
  if(is(mod,'glmmadmb')) {
    mod_fam <- mod$family
    mod_lnk <- mod$link
  }
  if(is(mod,'lmerMod')) {
    mod_fam <- family(mod)$family
    mod_lnk <- family(mod)$link
  }
  
  if(!exists('mod_fam')) stop(paste('Not implemented for',class(mod)[1],'objects!'))
  
  ### BINOMIAL FIT
  
  if (mod_fam=='binomial'|mod_fam=='binom') { 

    re_vars <- unlist(VarCorr(mod))
    fe_var <- var(model.matrix(mod) %*% fixef(mod))
    
    mod_check <- TRUE
  
    if(!exists('mod_lnk')) stop('Link function not defined? Not returned from model object.')
    
    if(mod_lnk=='logit') { 
      r2m <- fe_var/(fe_var+sum(re_vars)+(pi/3))
      r2c <- (fe_var+sum(re_vars))/(fe_var+sum(re_vars)+(pi/3))
    }
    
    if(mod_lnk=='probit') {
      r2m <- fe_var/(fe_var+sum(re_vars)+1)
      r2c <- (fe_var+sum(re_vars))/(fe_var+sum(re_vars)+1)
    }
    
    if(!exists('r2m')) stop(paste('Link function (',mod_lnk,') not defined!'))
  
    temp <- data.frame(Variance=round(c(fe_var, re_vars),4))
    row.names(temp) <- c('Fixed',names(re_vars))
  }
  
  ### POISSON FIT
  
  #   if (mod_fam=='poisson') { 
  #     
  #     mod_check <- TRUE
  #     
  #   }
  
  ### GAUSSIAN FIT
  
  if (mod_fam=='gaussian') { 
    
    re_vars <- unlist(VarCorr(mod))
    fe_var <- var(model.matrix(mod) %*% fixef(mod))
    
    r2m <- fe_var/(fe_var+sum(re_vars)+attr(VarCorr(mod),'sc')^2)
    r2c <- (fe_var+sum(re_vars))/(fe_var+sum(re_vars)+attr(VarCorr(mod),'sc')^2)
    
    mod_check <- TRUE
    
    temp <- data.frame(Variance=round(c(fe_var, re_vars, attr(VarCorr(mod),'sc')^2),4))
    row.names(temp) <- c('Fixed',as.data.frame(VarCorr(mod))$grp)
  }
  
  if(!exists('mod_check')) stop('Family not defined!')

  print(temp)
  #print(data.frame('Type'=c('Marginal (fixed)','Conditional (fixed+random)'),'R2'=round(c(r2m,r2c),4)),row.names=F)
  data.frame('R2'=round(c(r2m,r2c),4),row.names=c('Marginal','Conditional'))
}



# Plot semivariance (variogram) plot from variogram() data and add Loess smoother
# Bypasses horrible non-base graphic plot
#
# Argument ymax can set maximum value for y axis
plotvario <- function(vgdat, ymax=NULL, axes=TRUE) {
  if (is.null(ymax)) {
    ymax <- max(vgdat$gamma)
  }
  if(axes==TRUE) {
    plot(vgdat$dist, vgdat$gamma, ylim=c(0,ymax), 
         pch=16, xlab='Distance',ylab='Semivariance')
  }
  if(axes==FALSE) {
    plot(vgdat$dist, vgdat$gamma, ylim=c(0,ymax), 
         pch=16, xlab='Distance',ylab='Semivariance', yaxt='n', xaxt='n')
  }
  lines(vgdat$dist, loess(gamma~dist, data=vgdat)$fit, col='grey')
}

plotbubble <- function(x, y, z, size=1, xlab='x', ylab='y', alpha=NULL, axt='', colsign=TRUE) {
  if(colsign==TRUE) {
    cols <- factor(sign(z))
    levels(cols) <- c('grey','black')
    cols <- as.vector(cols)  
  } else {
    cols <- rep('grey', length(z))
  }
  if(axt=='n') {
    plot(x, y, col=cols, pch=16, cex=(abs(scale(z))+1)*size, xaxt='n',yaxt='n', xlab=xlab, ylab=ylab)
  } else {
    plot(x, y, col=cols, pch=16, cex=(abs(scale(z))+1)*size, xlab=xlab, ylab=ylab)    
  }
}

### Replicate arm:sim():

# This is a manual version of Gelman's sim() function.
# sim_man <- function(mod, S) {
#   sig.hat <- sigma(mod)   # Model estimated residual standard error
#   V.beta <- vcov(mod)     # Model estimated variance-covariance matrix
#   n <- nrow(as.data.frame(model.matrix(mod)))
#   k <- attr(logLik(mod),'df')
#   sigma <- as.vector(NULL)
#   beta <- as.data.frame(NULL)
#   for (s in 1:S) {
#     sigma <- c(sigma, sig.hat*sqrt((n-k)/rchisq(1,n-k)))
#     beta <- rbind(beta, mvrnorm(1, fixef(mod), V.beta*sigma[s]^2))
#   }
#   names(beta) <- names(fixef(mod))
#   return(list(fixef=beta, sigma=sigma))
# }
# 
# testmod <- glmer(OCC_PIPS ~ fSECTION*TURB + (1|SITE/TRSCT), data=bats_nona, 
#                  family=binomial(link='cloglog'))
# testmod <- standardize(testmod)
# 
# par(mfrow=c(3,4))
# par(mar=c(1,1,1,1))
# arm_sim <- density(attr(sim(testmod, 1000),'fixef')[,i])
# man_sim <- density(sim_man(testmod, 1000)$fixef[,i])
# for(i in 1:length(fixef(testmod))) {
#   ylims <- max(c(arm_sim$y, man_sim$y))
#   plot(arm_sim, main=names(fixef(testmod))[i])
#   lines(man_sim,col='red')
# }


### BEN BOLKER PREDICTION FUNCTION STUFF:

### http://rpubs.com/bbolker/glmmchapter

### Getting predicted values from an lme4 model (or an MCMCglmm model) is fairly straightforward: 
###  in this case by specifying re.form=NA we're saying that we want the population-level 
###  prediction, i.e. setting the random effects to zero and getting a prediction for an 
###  average (or unknown) block:

# pframe <- data.frame(ttt=factor(levels(culcita_dat$ttt),
# levels=levels(culcita_dat$ttt)))
# cpred1 <- predict(cmod_lme4_L,re.form=NA,newdata=pframe,type="response")

### Computing confidence intervals on the predicted values is relatively easy if we're 
### willing to completely ignore the random effects, and the uncertainty of the random effects. 
### Here is a generic function that extracts the relevant bits from a fitted model and returns 
###  the confidence intervals for predictions:
  
#   easyPredCI <- function(model,newdata,alpha=0.05) {
#     ## baseline prediction, on the linear predictor (logit) scale:
#     pred0 <- predict(model,re.form=NA,newdata=newdata)
#     ## fixed-effects model matrix for new data
#     X <- model.matrix(formula(model,fixed.only=TRUE)[-2],
#                       newdata)
#     beta <- fixef(model) ## fixed-effects coefficients
#     V <- vcov(model)     ## variance-covariance matrix of beta
#     pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
#     ## inverse-link (logistic) function: could also use plogis()
#     linkinv <- model@resp$family$linkinv
#     ## construct 95% Normal CIs on the link scale and
#     ##  transform back to the response (probability) scale:
#     crit <- -qnorm(alpha/2)
#     linkinv(cbind(lwr=pred0-crit*pred.se,
#                   upr=pred0+crit*pred.se))
#   }
# cpred1.CI <- easyPredCI(cmod_lme4_L,pframe)

