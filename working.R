mydat <- read.csv('boxplottest.csv')

boxPplot <- function(x, groups, scaler=0) {
  df <- data.frame(x=x, groups=groups)
  df$groups <- factor(df$groups)
  
  if(scaler==0) scaler <- nlevels(df$groups)*5  
  
  plot(1:nlevels(df$groups), type='n', 
       ylim=c(min(df$x),max(df$x)), 
       xlim=c(0.5,nlevels(df$groups)+0.5), 
       xlab='Group', ylab='y',
       xaxt='n')
  axis(1,at=1:nlevels(df$groups), labels=levels(df$groups))
  
  for(i in 1:nlevels(df$groups)) {
    y <- df$x[df$groups==levels(df$groups)[i]]
    freq <- hist(y, plot=F, breaks=length(y)*2*scaler)
    nbpb <- freq$counts[findInterval(y, freq$breaks)]-1
    jit_x <- rnorm(length(nbpb),0,nbpb/scale)
    points(jit_x+i, y)
  }
  
}