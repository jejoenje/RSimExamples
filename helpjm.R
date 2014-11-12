#helpjm.r

# JM HELPER FUNCTIONS

# Function to plot density curve of given vector dvect,
#  adding a vertical line at some reference point pval and naming the
#  plot 'pname'.
dplot <- function(dvect, pval, pname='') {
  d <- density(dvect)
  plot(d, main=pname)
  dmax <- max(d$y)
  lines(c(pval,pval),c(0,dmax),col='red',lty='dotted')
}

# Calculate overdispersion estimate a la Zuur 2009 (pp 224), from glm, glmer fit.
# This uses the ratio between the sum of the squared Pearson residuals (which is 
#   equivalent to the Pearson Chi-squared goodness of fit statistic) and the residual 
#   degrees of freedom, and should ideally be near 1.
# Note that this is a rough estimate only as arguably this underestimates the number
#   of "parameters" inferred by the random structure of any model.
dispZuur <- function(mod) {
  return(sum(resid(mod, type='pearson')^2)/(nrow(model.matrix(mod))-attr(logLik(mod),'df')))
}

# Plot semivariance (variogram) plot from variogram() data and add Loess smoother
# Bypasses horrible non-base graphic plot
#
# Argument ymax can set maximum value for y axis
plotvario <- function(vgdat, ymax=NULL) {
  if (is.null(ymax)) {
    ymax <- max(vgdat$gamma)
  }
  plot(vgdat$dist, vgdat$gamma, ylim=c(0,ymax), 
       pch=16, xlab='Distance',ylab='Semivariance')
  lines(vgdat$dist, loess(gamma~dist, data=vgdat)$fit, col='grey')
}

plotbubble <- function(x, y, z, size=1, xlab='x', ylab='y', alpha=NULL, axt='') {
  cols <- factor(sign(z))
  levels(cols) <- c('grey','black')
  cols <- as.vector(cols)
  if(axt=='n') {
    plot(x, y, col=cols, pch=16, cex=(abs(scale(z))+1)*size, xaxt='n',yaxt='n', xlab=xlab, ylab=ylab)
  } else {
    plot(x, y, col=cols, pch=16, cex=(abs(scale(z))+1)*size, xlab=xlab, ylab=ylab)    
  }
}

