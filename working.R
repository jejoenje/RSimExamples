mydat <- read.csv('boxplottest.csv')

boxPplot <- function(x, groups, scale=10) {
  df <- data.frame(x=x, groups=groups)
  df$groups <- factor(df$groups)
  
  for(i in 1:nlevels(df$groups)) {
    y <- df$x[df$groups==levels(df$groups)[i]]
    freq <- hist(y, plot=F, breaks=length(y)*10)
    nbpb <- freq$counts[findInterval(y, freq$breaks)]-1
    jit_x <- rnorm(length(nbpb),0,nbpb/scale)
  }
}