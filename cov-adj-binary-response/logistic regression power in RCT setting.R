# functions for power logistic regression

fun.d<-function(nsample, drug.allocation, 
                alpha,  beta.drug,
                seed=NULL){ 
  
  if (!is.null(seed)) set.seed(seed)
  
  drug<- (rbinom(nsample, 1, prob =drug.allocation ))   
  
  Xmat <- model.matrix(~ drug )
  beta.vec <- c(alpha,  beta.drug )
  
  lin.pred <- Xmat[,] %*% beta.vec                 # Value of lin.predictor
  exp.p <- exp(lin.pred) / (1 + exp(lin.pred))     # Expected proportion
  y <- rbinom(n = nsample, size = 1, prob = exp.p) # Add binomial noise
  #y<- runif(nsample) <  exp.p                     # alternatively ads noise in this way
  
  d<-as.data.frame(cbind(y, drug))         # create a dataset
  
  return(d)
  
}

# lrtest
simfunc <- function(d) {
  fit1 <- glm(y  ~ drug , d, family = binomial) 
  fit2 <- glm(y  ~ 1, d, family = binomial) 
  c(anova(fit1,fit2, test='Chisq')[2,5] )
}

# wald test

simfunc <- function(d) {
  fit1 <- glm(y  ~ drug , d, family = binomial) 
  c( summary(fit1)$coef["drug","Pr(>|z|)"] )
}


# (foo<-fun.d(nsample=20, drug.allocation=0.5,  
#             alpha=log(.3),  beta.drug=log(1)
#             , seed=NULL))

p1<-.35
odds <- .35/.65
or <- 1.5
prior.p<-p1
prior.odds<-prior.p/(1-prior.p)
post.odds<-prior.odds*or
p2<-post.p<-(post.odds/(1+post.odds))



out <- replicate(1000, simfunc(fun.d( nsample=300, drug.allocation=0.5,  
                                       alpha=log(odds),  
                                       beta.drug=log(or))))
mean( out < 0.05, na.rm=TRUE )        ##may get NA

bpower(p1=.35, odds.ratio=c(1.5 ), n=300, alpha=c(.05))


### just want to see what se is 


simfunc <- function(d) {
  fit1 <- glm(y  ~ drug , d, family = binomial) 
  c( summary(fit1)$coef["drug","Std. Error"] )
}



bsamsize(p1=.35, p2=p2, fraction=.5, alpha=.05, power=.9)

out <- replicate(1000, simfunc(fun.d( nsample=536*2, drug.allocation=0.5,  
                                      alpha=log(odds),  
                                      beta.drug=log(or))))

se. <- mean(out)



d <- fun.d( nsample=536*2, drug.allocation=0.5,  
            alpha=log(p1),  
            beta.drug=log(or))
table(d)/nrow(d)

#se
sqrt(sum(1/as.vector(table(d))))
sqrt(sum(1/as.vector(table(d)[,1])))
fit1 <- glm(y  ~ drug , d, family = binomial) 

summary(fit1)










#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# check this also
# https://www.dartmouth.edu/~eugened/
# We emphasize that the Wald test should be used to match a typically used coefficient significance testing. 
# https://www.dartmouth.edu/~eugened/
# frank Harrell
library(Hmisc)
bpower(p1=.35, odds.ratio=c(1.5,2), n=300, alpha=c(.05))
bpower.sim(p1=.35, odds.ratio=c(1.5), n=300, alpha=c(.05))


### just showing how to get  se 
tmp <- function (p1, p2, odds.ratio, percent.reduction, n, n1, n2, alpha = 0.05) 
{
  if (!missing(odds.ratio)) 
    p2 <- p1 * odds.ratio/(1 - p1 + p1 * odds.ratio)
  else if (!missing(percent.reduction)) 
    p2 <- p1 * (1 - percent.reduction/100)
  if (!missing(n)) {
    n1 <- n2 <- n/2
  }
  z <- qnorm(1 - alpha/2)
  q1 <- 1 - p1
  q2 <- 1 - p2
  pm <- (n1 * p1 + n2 * p2)/(n1 + n2)
  ds <- z * sqrt((1/n1 + 1/n2) * pm * (1 - pm))
  ex <- abs(p1 - p2)
  sd <- sqrt(p1 * q1/n1 + p2 * q2/n2)
  c(Power = 1 - pnorm((ds - ex)/sd) + pnorm((-ds - ex)/sd), se=sd)
}


p1<-.35
OR=1.5
p2 <- (.35/(1-.35)*1.5)
p2 <- p2/(1+p2)

Po <- bsamsize(p1,  p2, fraction=.5, alpha=.05, power=.9)

tmp(p1=.35, odds.ratio=c(1.5), n=535*2, alpha=c(.05))


sqrt(p1*(1-p1) / 535 +  p2*(1-p2)/535)


