### 22/11/2012. Works OK. Recovers parameters well.

MyWinBugsDir <- 'c:/WinBUGS14'

library(pscl)
library(MASS)
library(R2WinBUGS)
library(R2jags)

n.sims <- 10000

n.obs <- 1000
a <- 2
b <- 3
overdisp <- 0.75
psi <- 0.25         # Probability of a 'false' zero

x <- runif(n=n.obs, 0, 1)

dat <- data.frame(x=x); rm(x)

ZINB_est_a <- as.vector(NULL)
ZINB_est_b <- as.vector(NULL)
ZINB_est_psi <- as.vector(NULL)
ZINB_est_overdisp <- as.vector(NULL)

system.time(
{
  pb <- txtProgressBar(min = 0, max = n.sims, style = 3, width=getOption('width')/2)
  for(i in 1:n.sims) {
    W <- rbinom(n = n.obs, size=1, prob=1-psi)
    mu <- exp(a+b*dat$x)
    mu.eff <- W*mu

    dat$y <- rnbinom(n=n.obs, size=overdisp, prob=overdisp/(overdisp+mu.eff))   # For prob= see rnbinom help file, and Zuur et al 2009, p. 199.

    zinb.fit <- zeroinfl(y ~ x | 1, data=dat, dist='negbin')

    ZINB_est_a <- c(ZINB_est_a, as.numeric(coef(zinb.fit)['count_(Intercept)']))
    ZINB_est_b <- c(ZINB_est_b, as.numeric(coef(zinb.fit)['count_x']))
    ZINB_est_psi <- c(ZINB_est_psi, as.numeric(coef(zinb.fit)['zero_(Intercept)']))
    ZINB_est_overdisp <- c(ZINB_est_overdisp, zinb.fit$theta)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  par(mfrow=c(2,2))
  hist(ZINB_est_a)
  hist(ZINB_est_b)
  hist(ZINB_est_overdisp)
  hist(plogis(ZINB_est_psi))

  }

)


sink('NB_glm_test.txt')
cat('
  model {

    size ~ dgamma(0.001, 0.001)
    for (i in 1:2) { b[i] ~ dnorm(0.0, 0.001) }
    
    for (i in 1:N) {
      y[i] ~ dnegbin(p[i], size)
      p[i] <- size / (size + mu[i])
      log(mu[i]) <- max(-100, min(100, eta.mu[i]))
      eta.mu[i] <- b[1] + b[2]*x[i]
    }
  }  

    ', fill=T)
sink()

win.data <- list(
  N = length(dat$y),
  y = dat$y,
  x = dat$x
  )

inits <- function () {
  list(b           = rnorm(2, 0, 0.01),
       size        = 1
  )}

params <- c('b','size')

ni <- 5000
nc <- 3
nb <- 4000
nt <- 1


system.time(
  {
  #Start WINBUGS sampler
  out <- bugs(data = win.data,
              inits = inits,
              parameters = params,
              model = "NB_glm_test.txt",
              n.thin = nt,
              n.chains = nc,
              n.burnin = nb,
              n.iter = ni,
              debug = F,
              bugs.directory = MyWinBugsDir)
  
  par(mfrow=c(2,2))
  hist(out$sims.list$b[,1])
  hist(out$sims.list$b[,2])
  hist(out$sims.list$size)
  betas <- apply(out$sims.list$b,2,median)
  X <- X[order(X[,'x']),]
  pred_y <- exp(X %*% betas)
  plot(dat$x, dat$y)
  points(X[,'x'], pred_y, type='l')  
  }
  )

system.time (
{
  out_jags <- jags(
    data = win.data,
    parameters = params,
    inits = inits,
    model.file = "NB_glm_test.txt",
    n.thin = nt,
    n.chains = nc,
    n.burnin = nb,
    n.iter = ni
  )
}
  )


