#helpjm.r

# JM HELPER FUNCTIONS

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


