### Another simulated zero-inflated poisson example, from MY OWN BRAIN.
### JM 10/01/2013
###
### Wanted to clarify -
### (1) What the binomial part models. As shown below, the probability produced is
###     the probability that 'real' non-zero counts (from the poisson part) are 
###     'detected' (ie multiplied with a 1 rather than a 0).
### (2) The poisson model does indeed also predict zero's: these are 'real' zeros.


library(pscl)

### Number of fake obs:
runs <- 1000

N <- 100

### Probability psi = detection:

p_detect <- 0.25

### 1-psi = given actual presence (count > 0 below), species NOT detected:

p_missed <- 1-p_detect

### a is the average count
### b is the increase of count with x
a <- 1
b <- 2

ests <- as.data.frame(NULL)
for (i in 1:runs) {
  ### 100 fake observations:
  
  obs <- rbinom(N, prob=p_missed, size=1)
  
  ### 100 fake x's:
  x <- runif(N, 0, 1)
  
  mu <- exp(a+b*x)
  mu.eff <- mu*obs
  
  Y <- rpois(N, mu.eff)
  
  cbind(obs, Y)

  mod <- zeroinfl(Y ~ x | 1, dist='poisson')
  
  temp <- data.frame(a=coef(mod)[1], 
                     b=coef(mod)[2], 
                     p_detect=plogis(coef(mod)['zero_(Intercept)']),row.names=NULL)
  ests <- rbind(ests, temp)
}
apply(ests, 2, mean)
