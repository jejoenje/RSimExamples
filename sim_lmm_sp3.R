# sim_lmm_sp()
#
# Simulates normally distributed data (y or y_sp) as a linear function of continuous predictor x 
# (a=fixed intercept, b = fixed slope), and grouped by n_site "sites", each with n_obs_site 
# observations. In each site, normally distributed site-specific errors are added, as well as 
# residual error e. 
# Simulated response y is not spatially autocorrelated, whereas simulated response 
# y_sp is: for each observation, random coordinates xc (between xcmin and xcmax) and yc (between 
# ycmax and ycmin) are added, with the correlation strength modelled an exponential funciton of 
# distance and parameter rho. 
#
# sp = FALSE; no spatially correlated errors 
# spcor = TRUE; if T and sp = T, spatially correlated errors are added, which are correlated to the RE.
#  if F, spatially correlated errors are independent of RE.

sim_lmm_sp3 <- function(xcmin,           # Minimum x coordinate
                       xcmax,           # Maximum x coordinate
                       ycmin,           # Minimum y coordinate
                       ycmax,           # Maximum y coordinate
                       rho,             # Spatial correlation coefficient
                       re_sd,           # Random effect SD for SITE
                       e_sp_sd,         # Spatial error SD
                       e,               # Random (residual) error
                       a,               # Fixed intercept
                       b,               # Fixed slope
                       n_site,          # Number of sites ("groups")
                       n_obs_site,      # Number of observations per site
                       sp=FALSE,        # Spatially correlated errors?
                       spcor=FALSE,     # Spatially correlated errors clustered by RE? (Only when sp=T)
                       site_scale=100)  # Scale of between vs. within RE group SAC (Only used when sp=T)
{
  # Total number of observation is number of sites times number of obs per site:
  N <- n_site*n_obs_site
    
  if(sp==FALSE) {
  
    # First columns in data frame are random xc and yc coordinates for each obs, as well as a 
    #  measured continuous predictor x:
    mydat <- data.frame(xc=runif(N, xcmin, xcmax), yc=runif(N, ycmin, ycmax), x=runif(N, 0,5))
    
    # Calculate distance matrix between all points:
    d <- as.matrix(dist(cbind(mydat$xc, mydat$yc)))
    # m is the correlation matrix between all points, as function of distance and parameter rho:
    m <- exp(-rho*d)
    # e_sp are spatially correlated error terms for each observation, as draws from multivariate
    #  normal distribution with 0 means and covariance matrix m * error variance:
    e_sp <- mvrnorm(1, rep(0, N), m*e_sp_sd^2)
    
    # Now add some site errors. 
    # First add site identifiers to data
    mydat$site <- factor(rep(paste('site',1:n_site, sep=''),each=n_obs_site))
    # In this case, site-specific errors are entirely independent of the spatial correlation:
    site_e <- rnorm(n_site, 0, re_sd)       # Error in each site is N(0, re_sd)
    site_e <- rep(site_e, each=n_obs_site)  # Repeat these errors n_obs_site times for adding to simmed data.
    
    # Linear predictor is fixed + re + error.
    
    lin_y <- rnorm(N, a + b*mydat$x + site_e, e)
    
  } else {
    
    if(spcor==FALSE) {
      
      # Add spatially correlated errors, INdependent of RE:
      
      # First columns in data frame are random xc and yc coordinates for each obs, as well as a 
      #  measured continuous predictor x:
      mydat <- data.frame(xc=runif(N, xcmin, xcmax), yc=runif(N, ycmin, ycmax), x=runif(N, 0,5))
      
      # Calculate distance matrix between all points:
      d <- as.matrix(dist(cbind(mydat$xc, mydat$yc)))
      # m is the correlation matrix between all points, as function of distance and parameter rho:
      m <- exp(-rho*d)
      # e_sp are spatially correlated error terms for each observation, as draws from multivariate
      #  normal distribution with 0 means and covariance matrix m * error variance:
      e_sp <- mvrnorm(1, rep(0, N), m*e_sp_sd^2)
      
      # Now add some site errors. 
      # First add site identifiers to data
      mydat$site <- factor(rep(paste('site',1:n_site, sep=''),each=n_obs_site))
      # In this case, site-specific errors are entirely independent of the spatial correlation:
      site_e <- rnorm(n_site, 0, re_sd)       # Error in each site is N(0, re_sd)
      site_e <- rep(site_e, each=n_obs_site)  # Repeat these errors n_obs_site times for adding to simmed data.
      
      # Linear predictor is fixed + re + spatial error + error.
      
      lin_y <- rnorm(N, a + b*mydat$x + site_e + e_sp, e) 
      
    } else {
      
      # Add spatially correlated errors, dependent on RE:
      
      # Site "center" coordinates:
      sitex <- runif(n_site,xcmin,xcmax)
      sitey <- runif(n_site,ycmin,ycmax)
      # Individual obs coordinates are deviations from the site coordinates:
      x <- as.vector(NULL)
      y <- as.vector(NULL)
      for(i in 1:n_site) {
        x <- c(x,sitex[i]+rnorm(n_obs_site, 0, (xcmax-xcmin)/site_scale))
        y <- c(y,sitey[i]+rnorm(n_obs_site, 0, (ycmax-ycmin)/site_scale))
      }
      # First columns in data frame are random xc and yc coordinates for each obs, as well as a 
      #  measured continuous predictor x:
      mydat <- data.frame(xc=x, yc=y, x=runif(N, 0,5))
      
      # Calculate distance matrix between all points:
      d <- as.matrix(dist(cbind(mydat$xc, mydat$yc)))
      # m is the correlation matrix between all points, as function of distance and parameter rho:
      m <- exp(-rho*d)
      # e_sp are spatially correlated error terms for each observation, as draws from multivariate
      #  normal distribution with 0 means and covariance matrix m * error variance:
      e_sp <- mvrnorm(1, rep(0, N), m*e_sp_sd^2)
      
      # Now add some site errors. 
      # First add site identifiers to data
      mydat$site <- factor(rep(paste('site',1:n_site, sep=''),each=n_obs_site))
      # In this case, site-specific errors are entirely independent of the spatial correlation:
      site_e <- rnorm(n_site, 0, re_sd)       # Error in each site is N(0, re_sd)
      site_e <- rep(site_e, each=n_obs_site)  # Repeat these errors n_obs_site times for adding to simmed data.
      
      # Linear predictor is fixed + re + spatial error + error.      
      
      lin_y <- rnorm(N, a + b*mydat$x + site_e + e_sp, e) 
      
    }
  }
  

  # Observed y is some backtransformation (if a non-ID link is used) of lin_y's:
  mydat$y <- lin_y
  
  return(mydat)
}