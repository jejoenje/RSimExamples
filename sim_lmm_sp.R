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
# IMPORTANT NOTE: HERE, THE SPATIAL CORRELATION IS ENTIRELY INDEPEDENT OF THE RANDOM EFFECT OF SITE
# 

sim_lmm_sp <- function(xcmin,           # Minimum x coordinate
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
                       n_obs_site)      # Number of observations per site
{
  # Total number of observation is number of sites times number of obs per site:
  N <- n_site*n_obs_site
  
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
  
  # Linear predictors 
  lin_y <- a + b*mydat$x + site_e + rnorm(N, 0, e)           # Spatially correlated data
  lin_y_sp <- a + b*mydat$x + e_sp + site_e + rnorm(N, 0, e) # Spatially correlated data
  
  # Observed y is some backtransformation (if a non-ID link is used) of lin_y's:
  mydat$y <- lin_y
  mydat$y_sp <- lin_y_sp
  
  return(mydat)
}