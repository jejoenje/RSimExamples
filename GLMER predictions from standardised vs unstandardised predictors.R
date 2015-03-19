library(lme4)
library(arm)

test_us <- glmer(OCC_PIPS ~ SECTION + TURB + WINDS + (1|SITE/TRSCT), 
                 data=bats_nona, family=binomial)

pred_us <- data.frame(SECTION=1:5, 
                      TURB=factor(rep('single',5),levels=c('single','multiple')), 
                      WINDS=rep(mean(bats_nona$WINDS),5))

predict(test_us, newdata=pred_us, re.form=NA, type='response')

pred_us_m <- cbind(rep(1, 5), pred_us)
pred_us_m$TURB <- 0
pred_us_m <- as.matrix(pred_us_m)
plogis(pred_us_m %*% fixef(test_us))

# Plot predictions for SECTION:
barplot(predict(test_us, newdata=pred_us, re.form=NA, type='response'))

# Plot predictions for WIND:
windseq <- seq(min(bats_nona$WINDS), max(bats_nona$WINDS), 0.1)
pred_us2 <- data.frame(SECTION=rep(1, length(windseq)),
                       TURB=factor(rep('single',length(windseq)),levels=c('single','multiple')),
                       WINDS=windseq
                       )
predict(test_us, newdata=pred_us2, re.form=NA, type='response')
plot(windseq, predict(test_us, newdata=pred_us2, re.form=NA, type='response'), type='l')


### Attempt to retrieve the same point predictions from standardized model:

test_z <- standardize(test_us)
pred_z <- data.frame(z.SECTION=unique((bats_nona$SECTION-mean(bats_nona$SECTION))/(2*sd(bats_nona$SECTION))), 
                     c.TURB=rep(1-mean(as.numeric(bats_nona$TURB)), 5), 
                     z.WINDS=rep(0,5))
predict(test_z, newdata=pred_z, re.form=NA, type='response')
# Plot predictions for SECTION:
barplot(predict(test_z, newdata=pred_z, re.form=NA, type='response'))

windseq_z <- seq(min(rescale(bats_nona$WINDS)), max(rescale(bats_nona$WINDS)), 0.1)
pred_z2 <- data.frame(z.SECTION=rep((1-mean(bats_nona$SECTION))/(2*sd(bats_nona$SECTION)), length(windseq_z)),
                      c.TURB=rep(1-mean(as.numeric(bats_nona$TURB)), length(windseq_z)),
                      z.WINDS=windseq_z
                      )
predict(test_z, newdata=pred_z2, re.form=NA, type='response')
plot(windseq_z, predict(test_z, newdata=pred_z2, re.form=NA, type='response'), type='l')
# And with backtransformed axis:
plot(windseq_z*2*sd(bats_nona$WINDS)+mean(bats_nona$WINDS), 
     predict(test_z, newdata=pred_z2, re.form=NA, type='response'), type='l')

# So backtransformation is (y is standardised var and x raw var):
# x = y*2*sd(x)+mean(x)



