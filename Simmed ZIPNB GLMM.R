### 22/11/2012
### Works ok, but seems to massively underestimate the scale/overdispersion parameter for the neg. bin, 
###  as well as overestimate (?) the gamma1...


source('J:/035 Solway Firth/Data analysis/Operation Year 2 analysis/Latest stuff from Alain 081112/MCMCSupportHighstat.R')

MyWinBugsDir <- 'c:/WinBUGS14'

library(pscl)
library(MASS)
library(R2WinBUGS)
library(R2jags)
library(lattice)
library(lme4)


n.groups <- 50
n.obs.per.group <- 50
n.obs <- n.groups*n.obs.per.group
b1 <- 0
b2 <- 3
gamma1 <- 0.75
sigma.group <- 0.75
overdisp <- 1
a <- rnorm(n.groups, 0, sigma.group)

x <- runif(n=n.obs, 0, 1)
group_a <- rep(a, each=n.obs.per.group)
group <- gl(n.groups, n.obs.per.group)

dat <- data.frame(x=x, group_a=group_a, group=group); rm(x, group_a, group)

psi <- exp(gamma1)/(1+exp(gamma1))
dat$W <- rbinom(n.obs, size=1, prob=1-psi)

dat$lin.y <- b1+b2*dat$x+dat$group_a
dat$y <- dat$W*rnbinom(n=n.obs, 
                 size=overdisp, 
                 prob=(overdisp/(overdisp+exp(dat$lin.y))) )
                 # For prob= see rnbinom help file, and Zuur et al 2009, p. 199.
dat$y[is.nan(dat$y)] <- 0
xyplot(y~x|group, data=dat)
       
sink('ZIPNB_glmm_test.txt')
cat('
  model {

    size ~ dgamma(0.001, 0.001)
    tau.group <- 1/(sigma.group * sigma.group)
    sigma.group ~ dunif(0,10)
    for (i in 1:2) { b[i] ~ dnorm(0.0, 0.001) }
    for (i in 1:ngroups) { a[i] ~ dnorm(0.0, tau.group) }
    gamma1 ~ dnorm(0, 0.01)
    
    for (i in 1:N) {
      W[i] ~ dbern(psi.min)
      y[i] ~ dnegbin(mu.eff[i], size)
      mu.eff[i] <- W[i]*p[i]
      p[i] <- size / (size + mu[i])
      log(mu[i]) <- b[1] + b[2]*x[i] + a[group[i]]
    }
    psi.min <- 1-psi
    logit(psi) <- gamma1
  }  
    ', fill=T)
sink()

W <- dat$y
W[dat$y>0] <- 1

win.data <- list(
    N = length(dat$y),
    ngroups = nlevels(dat$group),
    y = dat$y,
    x = dat$x,
    group = as.numeric(dat$group)
  )

inits <- function () {
  list(a           = rnorm(nlevels(dat$group), 0, 0.01),
       b           = rnorm(2, 0, 0.01),
       size        = 1,
       sigma.group = 1,
       gamma1      = 1,
       W           = W
       
  )}

params <- c('a','b','gamma1','sigma.group','size')

ni <- 1000
nc <- 3
nb <- 500
nt <- 2


out <- bugs(data = win.data,
            inits = inits,
            parameters = params,
            model = "ZIPNB_glmm_test.txt",
            n.thin = nt,
            n.chains = nc,
            n.burnin = nb,
            n.iter = ni,
            debug = T,
            bugs.directory = MyWinBugsDir)

par(mfrow=c(2,2))
hist(out$sims.list$b[,1])
hist(out$sims.list$b[,2])
hist(out$sims.list$size)
hist(out$sims.list$sigma.group)

vars <- c('b[1]','b[2]','gamma1','sigma.group','size')
MyIterations(out, vars)
MyACF(out, vars)



