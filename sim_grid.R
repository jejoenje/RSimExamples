library(gstat)
library(raster)

### Function to simulate a 'landscape grid' on which n_obs observations fall.
### Use of gstat() based on http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/
### Generates square x-y landscape with sides equal to sqrt(n_obs+1), to make sure
###  all observations would fit onto unique cells if necessary.
###
### Uses gstat() to simulate unconditional random fields using an exponential VGM model.
sim_grid <- function(n_obs, sill=0.1, range=5, plotme=F) {
  
  xy <- expand.grid(xc=1:ceiling(sqrt(n_obs+1)), yc=1:ceiling(sqrt(n_obs+1)))
  spmod <- gstat(formula=z~1+xc+yc, locations=~xc+yc, dummy=T, 
                 beta=c(0,-0.05,0.05), model=vgm(psill=sill,range=range,model='Exp'), nmax=20)
  sp <- predict(spmod, newdata=xy, nsim=1)
  if (plotme==T) {
    plot_grid(sp)
  }
  return(sp)
}

### Function to plot a dataframe with x,y coordinates and some value sim1 as a raster:
plot_grid <- function(sp) {
  rast_sp <- raster(ncol=max(sp$x), nrow=max(sp$y), ymn = min(sp$y), ymx = max(sp$y), 
                    xmn = min(sp$x), xmx = max(sp$x))
  res(rast_sp) <- 1
  rast_sp[cellFromXY(rast_sp, cbind(sp$x,sp$y))] <- sp$sim1
  plot(rast_sp, col=brewer.pal(11,"RdGy"))
}
