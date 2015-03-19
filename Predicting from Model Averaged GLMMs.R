### Predicting from model-averaged results

library(lme4)
library(MuMIn)

load('m2z.Rdata')
load('m2z_set1.Rdata')
load('m2z_av.Rdata')
load('m2z_set1_mods.Rdata')

### So, m2z is a Binomial GLMM with cloglog link, a nested RE and a bunch of FE's:
summary(m2z)

### m2z_set1 is the output from a model selection procedure on m2z with some restrictions.
subset(m2z_set1, delta<4)

### m2z_av is the output from a model averaging procedure on the top delta<4 models above:
summary(m2z_av)

### So how do we predict from these averaged models / average paramter estimates?

### This extracts the link function used: the same for all models in the set.
linkinv <- m2z@resp$family$linkinv

### This makes a vector of one of two levels of c.TURB to predict for:
### (Originally this was a factor with two levels - because of the standardization method used (standardize() in arm), this
### results in the values being centred only so that the difference between the two factor levels is 1:)
c.turb <- unique(m2z@frame$c.TURB)[order(unique(m2z@frame$c.TURB))]

### First, check we can make OK predictions from the full model with predict() as well as manually.
### First, using predict()
### ... for the first value in c.turb:
pred_full_single <- predict(m2z, type='response', re.form=NA, newdata=data.frame(
  fSECTION=factor(1:5),
  c.TURB=rep(c.turb[1],5),
  z.MINTEMP=rep(0,5),
  z.DAYNO=rep(0,5),
  z.TTMIDN=rep(0,5),
  z.WINDS=rep(0,5),
  z.EDGED=rep(0,5),
  z.pTREE=rep(0,5),
  AREA_ha=rep(1,5)
))
### ... for the second value in c.turb:
pred_full_multiple <- predict(m2z, type='response', re.form=NA, newdata=data.frame(
  fSECTION=factor(1:5),
  c.TURB=rep(c.turb[2],5),
  z.MINTEMP=rep(0,5),
  z.DAYNO=rep(0,5),
  z.TTMIDN=rep(0,5),
  z.WINDS=rep(0,5),
  z.EDGED=rep(0,5),
  z.pTREE=rep(0,5),
  AREA_ha=rep(1,5)
))

### Now repeat the above predictions manually:
### Note that the prediction 'frame' looks different as we need to specify values for each parameter.
### ... for the first value in c.turb:
p_full_single <- cbind(
  rep(1,5),           # Intercept
  c(0,1,0,0,0),       # fSECTION2
  c(0,0,1,0,0),       # fSECTION3
  c(0,0,0,1,0),       # fSECTION4
  c(0,0,0,0,1),       # fSECTION5
  rep(c.turb[1],5),   # c.TURB [SINGLE]
  rep(0, 5),          # z.MINTEMP 
  rep(0, 5),          # z.DAYNO
  rep(0, 5),          # z.TTMIDN
  rep(0, 5),          # z.TTMIDN^2
  rep(0, 5),          # z.WINDS
  rep(0, 5),          # z.EDGED
  rep(0, 5),          # z.pTREE
  cbind(
    c(0,1,0,0,0),       # fSECTION2*c.TURB
    c(0,0,1,0,0),       # fSECTION3*c.TURB
    c(0,0,0,1,0),       # fSECTION4*c.TURB
    c(0,0,0,0,1)        # fSECTION5*c.TURB
  )*c.turb[1]
)
pred_manual_single <- linkinv(p_full_single %*% fixef(m2z))

### ... for the second value in c.turb:
p_full_mult <- cbind(
  rep(1,5),           # Intercept
  c(0,1,0,0,0),       # fSECTION2
  c(0,0,1,0,0),       # fSECTION3
  c(0,0,0,1,0),       # fSECTION4
  c(0,0,0,0,1),       # fSECTION5
  rep(c.turb[2],5),   # c.TURB [SINGLE]
  rep(0, 5),          # z.MINTEMP 
  rep(0, 5),          # z.DAYNO
  rep(0, 5),          # z.TTMIDN
  rep(0, 5),          # z.TTMIDN^2
  rep(0, 5),          # z.WINDS
  rep(0, 5),          # z.EDGED
  rep(0, 5),          # z.pTREE
  cbind(
    c(0,1,0,0,0),       # fSECTION2*c.TURB
    c(0,0,1,0,0),       # fSECTION3*c.TURB
    c(0,0,0,1,0),       # fSECTION4*c.TURB
    c(0,0,0,0,1)        # fSECTION5*c.TURB
  )*c.turb[2]
)
pred_manual_multiple <- linkinv(p_full_mult %*% fixef(m2z))

### Now check that these predictions are the same:
### First value in c.turb:
pred_full_single
pred_manual_single
pred_full_multiple
pred_manual_multiple
### These are identical as they should be.

### 
### So, we can use predict() on the full GLMM object.
### But how do we predict from model-averaged coefficients?
### First, let's see what happens if we use the predict() function on the averaging object:
### ... for the first value in c.turb:
pred_av_single <- predict(m2z_av, type='response', re.form=NA, newdata=data.frame(
  fSECTION=factor(1:5),
  c.TURB=rep(c.turb[1],5),
  z.MINTEMP=rep(0,5),
  z.DAYNO=rep(0,5),
  z.TTMIDN=rep(0,5),
  z.WINDS=rep(0,5),
  z.EDGED=rep(0,5),
  z.pTREE=rep(0,5),
  AREA_ha=rep(1,5)
))
### ... and for the second value:
pred_av_multiple <- predict(m2z_av, type='response', re.form=NA, newdata=data.frame(
  fSECTION=factor(1:5),
  c.TURB=rep(c.turb[2],5),
  z.MINTEMP=rep(0,5),
  z.DAYNO=rep(0,5),
  z.TTMIDN=rep(0,5),
  z.WINDS=rep(0,5),
  z.EDGED=rep(0,5),
  z.pTREE=rep(0,5),
  AREA_ha=rep(1,5)
))

### Now lets try to replicate the above results manually, by extracting the relevant parameters and
### multiplying with a prediction frame as above.
### However, which parameters should we use. There are two options available in the summary(m2z_av):
### 'Natural averaging', the first set of parameters in the summary(),
coefTable(m2z_av, full=FALSE)
### ... which basically average only over the models that include the parameter in question.
### Alternatively, there are the 'zero method' coefficients (the second set of coefficients in summary()).
### These correspond to the parameter estimates "with shrinkage" - they average all parameters over all
###  models, effectively implying that the parameter is zero in those models where the parameter isn't included.
### This has the effect of 'shrinking' the parameter estimate to zero:
coefTable(m2z_av, full=TRUE)

### So which should we use, and importantly, which does predict.averaging() use?
### Let's try to make manual predictions with each.
### Let's start with the natural averaged ones.
### ... for the first value in c.turb:
p_single <- cbind(
  rep(1,5),           # Intercept
  c(0,1,0,0,0),       # fSECTION2
  c(0,0,1,0,0),       # fSECTION3
  c(0,0,0,1,0),       # fSECTION4
  c(0,0,0,0,1),       # fSECTION5
  rep(0,5),           # z.DAYNO
  rep(0,5),           # z.EDGED
  rep(0,5),           # z.pTREE  
  rep(0,5),           # z.TTMIDN
  rep(0,5),           # z.TTMIDN^2
  rep(0,5),           # z.WINDS
  rep(c.turb[1],5),   # c.TURB
  rep(0,5),           # z.MINTEMP
  cbind(
    c(0,1,0,0,0),       # fSECTION2*c.TURB
    c(0,0,1,0,0),       # fSECTION3*c.TURB
    c(0,0,0,1,0),       # fSECTION4*c.TURB
    c(0,0,0,0,1)        # fSECTION5*c.TURB
  )*c.turb[1]
)
pred_natav_single <- linkinv(p_single %*% coefTable(m2z_av, full=FALSE)[,1])

### ... for the second value in c.turb:
p_multp <- cbind(
  rep(1,5),           # Intercept
  c(0,1,0,0,0),       # fSECTION2
  c(0,0,1,0,0),       # fSECTION3
  c(0,0,0,1,0),       # fSECTION4
  c(0,0,0,0,1),       # fSECTION5
  rep(0,5),           # z.DAYNO
  rep(0,5),           # z.EDGED
  rep(0,5),           # z.pTREE
  rep(0,5),           # z.TTMIDN
  rep(0,5),           # z.TTMIDN^2
  rep(0,5),           # z.WINDS
  rep(c.turb[2],5),   # c.TURB
  rep(0,5),           # z.MINTEMP
  cbind(
    c(0,1,0,0,0),       # fSECTION2*c.TURB
    c(0,0,1,0,0),       # fSECTION3*c.TURB
    c(0,0,0,1,0),       # fSECTION4*c.TURB
    c(0,0,0,0,1)        # fSECTION5*c.TURB
  )*c.turb[2]
)
pred_natav_multiple <- linkinv(p_multp %*% coefTable(m2z_av, full=FALSE)[,1])

### Now check these against the predictions from predict():
pred_av_single
pred_natav_single 
pred_av_multiple
pred_natav_multiple 
### So, predict.averaging() clearly does NOT use natural averaging parameters.

### What about the "shrunk", zero method, parameter estimates?
### The prediction frame is the same as used by the natural averaging calculation above, so:
### ... for the first value in c.turb:
(pred_zeroav_single <- linkinv(p_single %*% coefTable(m2z_av, full=TRUE)[,1]))
pred_av_single
(pred_zeroav_multiple <- linkinv(p_multp %*% coefTable(m2z_av, full=TRUE)[,1]))
pred_av_multiple
### While these point predictions are a lot closer than those from the natural averaged parameters,
###  its still not quite the same!

### So, the last option is that predict.averaging() actually calculates predictions from each individual
### model in the set individually, and then averages the predictions, rather than directly use averaged
### parameter estimates to make the predictions.

### So lets loop through each of the models in our candidate set (note that we need to have a dredge() object
###  fit with fit=T so that it contains the actual models rather than just the model selection output).
### In the following, I first make an output data frame, and then use the predict() function on each of the individual models
###  in the candidate set. For each, I store the five point predictions in the output frame.
### ... for the first value in c.turb:
allpreds_single <- as.data.frame(NULL)
for(i in 1:length(m2z_set1_mods)) {
  allpreds_single <- rbind(allpreds_single, predict(m2z_set1_mods[[i]], type='response', re.form=NA, newdata=data.frame(
    fSECTION=factor(1:5),
    c.TURB=rep(c.turb[1],5),
    z.MINTEMP=rep(0,5),
    z.DAYNO=rep(0,5),
    z.TTMIDN=rep(0,5),
    z.WINDS=rep(0,5),
    z.EDGED=rep(0,5),
    z.pTREE=rep(0,5),
    AREA_ha=rep(1,5)
  )))
}
### ... for the second value in c.turb:
allpreds_multp <- as.data.frame(NULL)
for(i in 1:length(m2z_set1_mods)) {
  allpreds_multp <- rbind(allpreds_multp, predict(m2z_set1_mods[[i]], type='response', re.form=NA, newdata=data.frame(
    fSECTION=factor(1:5),
    c.TURB=rep(c.turb[2],5),
    z.MINTEMP=rep(0,5),
    z.DAYNO=rep(0,5),
    z.TTMIDN=rep(0,5),
    z.WINDS=rep(0,5),
    z.EDGED=rep(0,5),
    z.pTREE=rep(0,5),
    AREA_ha=rep(1,5)
  )))
}

### Now make a little function that can easily calculate weighted means using our model weights:
wmean <- function(x, w=NULL) {
  if(is.null(w)) { w <- rep(1/length(x), length(x)) }
  return(as.vector(as.numeric(sum(x*w))))
}

### Extract the model weights from the model averaging object:
modweights <- subset(m2z_set1, delta<4)$weight

### Use this function to calculate the weighted means of each of the ten predictions for each of the five
###  values, for both sets of predictions, and compare the weighted point predictions to those from predict.averaging():
(pred_avmanual_single <- as.vector(apply(allpreds_single, 2, function(x) wmean(x, modweights))))
pred_av_single
(pred_avmanual_multiple <- as.vector(apply(allpreds_multp, 2, function(x) wmean(x, modweights))))
pred_av_multiple
### These are exactly the same - hurrah!

### So predict.averaging() indeed calculates predictions for each model in the set individually, and produces a weighted
###  average of predictions summing model prediction * model weight. Note that here I averaged the 
###  BACK-TRANSFORMED predictions instead of first averaging on the link-scale and then backtransforming. This must clearly be
###  the default for predict.averaging() - which it is as per the helpfile 'backTransform=FALSE'.

### As per the above, this clearly gives DIFFERENT RESULTS to all other options. Thus, when model averageing, we have a bunch of
###  potential options when calculating predictions:
### 1. Predict from "natural average" estimates
### 2. Predict from "zero method" estimates
### 3. Predict from all individual models, backtransform, and average overall models weighted by model weight.
### 4. Predict from the "best" model (effectively making the averaging redundant)
### 5. Predict from the "full" model (effectively making the whole IT approach redundant)

### Using the output from some of the above calculations, here is a plot of the predictions calculated
###  using approach 1, 2, 3, and 5 listed above:
singles <- cbind(pred_natav_single, pred_zeroav_single, pred_av_single, pred_full_single)
multiples <- cbind(pred_natav_multiple, pred_zeroav_multiple, pred_av_multiple, pred_full_multiple)
dimnames(singles)[[2]] <- c("Natural", "Zero", "Prediction", "None (full model)")
dimnames(multiples)[[2]] <- c("Natural", "Zero", "Prediction", "None (full model)")
par(mfrow=c(2,1))
barplot(singles, beside=T, xlab='Averaging method', col=rep('grey',5))
barplot(multiples, beside=T, xlab='Averaging method', col=rep('grey',5))

### Note that the predictions from the predict.averaging() - i.e. a weighted average of predictions from
###  individual models - and those from using the 'zero' method parameter estimates are very similar indeed.
### The fact that they are not exactly the same begs the question whether there are situations (e.g. strong interactions, 
###  unstable model selection, colinear variables) where this might not be the case... it would be very interesting
###  indeed to know when exactly this would be the case.
