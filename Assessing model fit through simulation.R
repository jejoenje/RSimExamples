setwd("~/Documents/docs/000_R/")
library(foreign)
library(MASS)
library(ggplot2)
data <- read.dta("repdata.dta")

data <- data[data$onset != 4,]

data$warl <- as.factor(data$warl)
data$gdpenl <- as.numeric(as.vector(data$gdpenl))
data$lpopl <- as.numeric(as.vector(data$lpopl))
data$ncontig <- as.factor(data$ncontig)
data$Oil <- as.factor(data$Oil)
data$nwstate <- as.factor(data$nwstate)
data$instab <- as.factor(data$instab)

#estimate the model
model <- glm(onset ~ warl + gdpenl + lpopl + lmtnest + 
             ncontig + Oil + nwstate + instab + polity2l + 
             ethfrac, 
             data = data, family = "binomial")
summary(model)

#create a matrix with the variables used in the model and remove 
# rows with missing values
X <- cbind(1, data$warl, data$gdpenl, data$lpopl, 
              data$lmtnest, data$ncontig, data$Oil, 
              data$nwstate, data$instab, data$polity2l, data$ethfrac)

X <- na.omit(X)

#get the beta coefficients and variance-covariance matrix from the model
beta <- coef(model)
covvar <- vcov(model)

#create a matrix to hold the fake data sets
n <- nrow(X)
out <- matrix(nr = n, nc = 1000)

#simulate a 1000 fake data sets
for (s in 1:1000){
  b <- mvrnorm(1,beta,covvar)
  xb <- X %*% b
  p <- 1/(1 + exp(-xb))
  y.fake <- rbinom(n,1,p)
  out[,s] <- y.fake
}

#for each fake data set calculate the proportion of 1s
pred <- apply(out,2,function(x){
  p <- sum(x == 1)/length(x)
  return(p)
})

#get the proportion from the observed data
p.true <- sum(data$onset == 1)/length(data$onset)

#put the simulated proportions in a data frame for ggplot2
pred <- data.frame(pred = pred)


#plot the density of the simulated proportions and mark the observed proportion
plot <- ggplot(pred, aes(x = pred)) +
  geom_density() +
  theme_bw() +
  geom_segment(aes(x = p.true, xend = p.true, y = 0, yend = 100), 
               colour = "red", 
               arrow=arrow(length=unit(0.3,"cm"), 
                           ends = "first"))

plot


custom_sims <- as.vector(NULL)
for (i in 1:1000) {
  simmed_vector <- rbinom(n, 1, plogis(X %*% mvrnorm(1,beta,covvar)))
  custom_sims[i] <- sum(simmed_vector == 1)/length(simmed_vector)
}

custom_sims2 <- as.vector(NULL)
for (i in 1:1000) {
  pred <- predict(model, type="response")
  simmed_vector <- rbinom(length(pred),1,pred)
  custom_sims2[i] <- sum(simmed_vector == 1)/length(simmed_vector)
}

r_sims <- as.vector(NULL)
for (i in 1:1000) {
  r_vector <- simulate(model)[,1]
  r_sims[i] <- sum(r_vector == 1)/length(r_vector)
}

plotdata <- data.frame(custom_sims, r_sims)
plotdata2 <- data.frame(custom_sims2, r_sims)


test <- ggplot(plotdata, aes(x=custom_sims, y=r_sims)) + geom_point()
test

test2 <- ggplot(plotdata2, aes(x=custom_sims2, y=r_sims, xlim=c(0.01,0.023), ylim=c(0.01,0.023))) + geom_point()
test2





