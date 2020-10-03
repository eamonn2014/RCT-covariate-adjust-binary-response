set.seed(1)
# https://stats.stackexchange.com/questions/35940/simulation-of-logistic-regression-power-analysis-designed-experiments
repetitions = 1000
N = 10000
n = N/8
var1  = c(   .03,    .03,    .03,    .03,    .06,    .06,    .09,   .09)
var2  = c(     0,      0,      0,      1,      0,      1,      0,     1)
rates = c(0.0025, 0.0025, 0.0025, 0.00395, 0.003, 0.0042, 0.0035, 0.002)

var1    = rep(var1, times=n)
var2    = rep(var2, times=n)
var12   = var1**2
var1x2  = var1 *var2
var12x2 = var12*var2

significant = matrix(nrow=repetitions, ncol=7)

startT = proc.time()[3]
for(i in 1:repetitions){
  responses          = rbinom(n=N, size=1, prob=rates)
  model              = glm(responses~var1+var2+var12+var1x2+var12x2, 
                           family=binomial(link="logit"))
  significant[i,1:5] = (summary(model)$coefficients[2:6,4]<.05)
  significant[i,6]   = sum(significant[i,1:5])
  modelDev           = model$null.deviance-model$deviance
  significant[i,7]   = (1-pchisq(modelDev, 5))<.05
}
endT = proc.time()[3]
endT-startT

sum(significant[,1])/repetitions 

#######################

mydat <- data.frame( v1 = rep( c(3,6,9), each=2 ),
                     v2 = rep( 0:1, 3 ), 
                     resp=c(0.0025, 0.00395, 0.003, 0.0042, 0.0035, 0.002) )

fit0 <- glm( resp ~ poly(v1, 2, raw=TRUE)*v2, data=mydat,
             weight=rep(100000,6), family=binomial)
b0 <- coef(fit0)


simfunc <- function( beta=b0, n=10000 ) {
  w <- sample(1:6, n, replace=TRUE, prob=c(3, rep(1,5)))
  mydat2 <- mydat[w, 1:2]
  eta <- with(mydat2,  cbind( 1, v1, 
                              v1^2, v2,
                              v1*v2,
                              v1^2*v2 ) %*% beta )
  p <- exp(eta)/(1+exp(eta))
  mydat2$resp <- rbinom(n, 1, p)
  
  fit1 <- glm( resp ~ poly(v1, 2)*v2, data=mydat2,
               family=binomial)
  fit2 <- update(fit1, .~ poly(v1,2) )
  anova(fit1,fit2, test='Chisq')[2,5]
}

out <- replicate(100, simfunc(b0, 10000))
mean( out <= 0.05 )
hist(out)
abline(v=0.05, col='lightgrey')