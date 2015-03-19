### 23/11/2012. pscl (zeroinfl) fit works very well, recovers parameters perfectly.
###             WinBUGS fit seems faulty, suspect this is not correctly parameterised.

MyWinBugsDir <- 'c:/WinBUGS14'

library(pscl)
library(MASS)
library(R2WinBUGS)
library(R2jags)

n.sims <- 100

n.obs <- 500
a <- 2
b <- 3
c <- -1       # Int in ZERO
d <- 5        # Slope in ZERO
overdisp <- 0.75

x <- runif(n=n.obs, 0, 1)

dat <- data.frame(x=x); rm(x)

ZINB_est_a <- as.vector(NULL)
ZINB_est_b <- as.vector(NULL)
ZINB_est_c <- as.vector(NULL)
ZINB_est_d <- as.vector(NULL)
ZINB_est_overdisp <- as.vector(NULL)

BUGS_est_a <- as.vector(NULL)
BUGS_est_b <- as.vector(NULL)
BUGS_est_c <- as.vector(NULL)
BUGS_est_d <- as.vector(NULL)
BUGS_est_overdisp <- as.vector(NULL)

med_BUGS_est_a <- as.vector(NULL)
med_BUGS_est_b <- as.vector(NULL)
med_BUGS_est_c <- as.vector(NULL)
med_BUGS_est_d <- as.vector(NULL)
med_BUGS_est_overdisp <- as.vector(NULL)


sink('ZINB_glm_test.txt')
cat('
  model {

    size ~ dgamma(0.001, 0.001)
    for (i in 1:2) { beta[i] ~ dnorm(0.0, 0.001) }
    for (i in 1:2) { gamma[i] ~ dnorm(0.0, 0.001) }
    
    for (i in 1:N) {
      W[i] ~ dbern(psi.min[i])
      psi.min[i] <- 1-psi[i]
      logit(psi[i]) <- eta.psi[i]
      eta.psi[i] <- gamma[1] + gamma[2]*x[i]

      y[i] ~ dnegbin(p[i], size)
      p[i] <- size / (size + mu.eff[i])
      mu.eff[i] <- W[i]*mu[i]
      log(mu[i]) <- eta.mu[i]
      eta.mu[i] <- beta[1] + beta[2]*x[i]
    }
  }  

    ', fill=T)
sink()

ni <- 3000
nc <- 3
nb <- 2700
nt <- 3

system.time(
{
  pb <- txtProgressBar(min = 0, max = n.sims, style = 3, width=getOption('width')/2)
  for(i in 1:n.sims) {
    eta.psi <- c + d*dat$x
    psi <- exp(eta.psi)/(1+exp(eta.psi))
    W <- rbinom(n = n.obs, size=1, prob=1-psi)
    mu <- exp(a+b*dat$x)
    mu.eff <- W*mu

    dat$y <- rnbinom(n=n.obs, size=overdisp, prob=overdisp/(overdisp+mu.eff))   # For prob= see rnbinom help file, and Zuur et al 2009, p. 199.

    zinb.fit <- zeroinfl(y ~ x | x, data=dat, dist='negbin')

    ZINB_est_a <- c(ZINB_est_a, as.numeric(coef(zinb.fit)['count_(Intercept)']))
    ZINB_est_b <- c(ZINB_est_b, as.numeric(coef(zinb.fit)['count_x']))
    ZINB_est_c <- c(ZINB_est_c, as.numeric(coef(zinb.fit)['zero_(Intercept)']))
    ZINB_est_d <- c(ZINB_est_d, as.numeric(coef(zinb.fit)['zero_x']))    
    ZINB_est_overdisp <- c(ZINB_est_overdisp, zinb.fit$theta)
    
    W <- dat$y
    W[dat$y > 0] <- 1
    
    win.data <- list(
      N = length(dat$y),
      y = dat$y,
      x = dat$x,
      W = W)
    
    inits <- function () {
      list(beta        = rnorm(2, 0, 0.01),
           gamma       = rnorm(2, 0, 0.01),
           size        = 1
      )}
    
    params <- c('beta','gamma','size')
    
    zipnbfit <- bugs(data = win.data,
                     inits = inits,
                     parameters = params,
                     model = "ZINB_glm_test.txt",
                     n.thin = nt,
                     n.chains = nc,
                     n.burnin = nb,
                     n.iter = ni,
                     debug = F,
                     bugs.directory = MyWinBugsDir)
    
    BUGS_est_a <- c(BUGS_est_a, zipnbfit$mean$beta[1])
    BUGS_est_b <- c(BUGS_est_b, zipnbfit$mean$beta[2])
    BUGS_est_c <- c(BUGS_est_c, zipnbfit$mean$gamma[1])
    BUGS_est_d <- c(BUGS_est_d, zipnbfit$mean$gamma[2])
    BUGS_est_overdisp <- c(BUGS_est_overdisp, zipnbfit$mean$size)
    
    med_BUGS_est_a <- c(med_BUGS_est_a, zipnbfit$median$beta[1])
    med_BUGS_est_b <- c(med_BUGS_est_b, zipnbfit$median$beta[2])
    med_BUGS_est_c <- c(med_BUGS_est_c, zipnbfit$median$gamma[1])
    med_BUGS_est_d <- c(med_BUGS_est_d, zipnbfit$median$gamma[2])
    med_BUGS_est_overdisp <- c(med_BUGS_est_overdisp, zipnbfit$median$size)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  par(mfrow=c(3,2))
  hist(BUGS_est_a)
  hist(BUGS_est_b)
  hist(BUGS_est_c)
  hist(BUGS_est_d)  
  hist(BUGS_est_overdisp)

  }

)








#Start WINBUGS sampler

par(mfrow=c(3,2))
hist(zipnbfit$sims.list$beta[,1])
hist(zipnbfit$sims.list$beta[,2])
hist(zipnbfit$sims.list$gamma[,1])
hist(zipnbfit$sims.list$gamma[,2])
hist(zipnbfit$sims.list$size)
  


# system.time (
# {
#   out_jags <- jags(
#     data = win.data,
#     parameters = params,
#     inits = inits,
#     model.file = "NB_glm_test.txt",
#     n.thin = nt,
#     n.chains = nc,
#     n.burnin = nb,
#     n.iter = ni
#   )
# }
  )


