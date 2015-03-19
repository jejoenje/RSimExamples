### Implementation of a zero-inflated Poisson (ZIP) in WinBUGS and JAGS, using simulated data

### Example from Zuur et al 2012, p. 57-60.

### 21/11/2012. Works OK. Recovers parameters well.

library(R2WinBUGS)

setwd('~/Documents/docs/000_R/fake-data trials')

set.seed(12345)
beta1 <- 2
beta2 <- 2.5
gamma1 <- 2
N <- 250
X <- runif(N, min=0, max=1)

psi <- exp(gamma1)/(1+exp(gamma1))
W <- rbinom(N, size=1, prob=1-psi)

table(W)

mu <- exp(beta1 + beta2 * X)
mu.eff <- W * mu
Y <- rpois(N, lambda=mu.eff)

T1 <- glm(Y ~ X , family='poisson')
summary(T1)

deviance(T1)/T1$df.res

library(pscl)

Z1 <- zeroinfl(Y ~ X | 1)
summary(Z1)

win.data <- list(Y = Y, X = X, N = length(Y))

sink('ziptest.txt')
cat('
    model {

  # Priors
    beta[1] ~ dnorm(0, 0.001)
    beta[2] ~ dnorm(0, 0.001)
    gamma1 ~ dnorm(0, 0.001)

  # Likelihood
    for (i in 1:N) {
    # Binary part
      W[i] ~ dbern(psi.min1)
    # Count process
      Y[i] ~ dpois(mu.eff[i])
      mu.eff[i] <- W[i] * mu[i]
      log(mu[i]) <- beta[1] + beta[2] * X[i]
    }
    psi.min1 <- 1-psi
    logit(psi) <- gamma1

    }
    
    ', fill=T)
sink()

W <- Y
W[Y>0] <- 1

inits <- function() {
  list(beta = rnorm(2),
       gamma1 = rnorm(1),
       W = W
       )
}

params <- c('beta','gamma1')

nc <- 3
ni <- 1000
nb <- 500
nt <- 2
WinBugsDir = 'c:/WinBUGS14'

system.time( ZipSim1 <- bugs(data = win.data, 
                inits = inits, 
                parameters = params, 
                model = 'ziptest.txt', 
                n.thin=nt, 
                n.chains=nc, 
                n.burnin = nb, 
                n.iter = ni, 
                debug=F,
                bugs.directory=WinBugsDir) )

system.time( ZipSim1 <- bugs(data = win.data, 
                             inits = inits, 
                             parameters = params, 
                             model = 'ziptest.txt', 
                             n.thin=nt, 
                             n.chains=nc, 
                             n.burnin = nb, 
                             n.iter = ni, 
                             debug=F)
             
system.time ( ZipSim1.jags <- jags(data = win.data,
                inits = inits,
                parameters = params,
                model.file = 'ziptest.txt',
                n.thin=nt, 
                n.chains=nc, 
                n.burnin = nb, 
                n.iter = ni
                ) )