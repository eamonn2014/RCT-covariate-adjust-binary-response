# R code to help reproduce app
# logistic regression covariate adjustment

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    rm(list=ls()) 
    set.seed(333) # reproducible
    
    library(Hmisc)
    library(reshape)
    library(rms)
    library(ggplot2)
    library(tidyverse)
    library(Matrix)
    
    options(max.print=1000000)    
    
    
    ## convenience functions
    p0 <- function(x) {formatC(x, format="f", digits=0)}
    p1 <- function(x) {formatC(x, format="f", digits=1)}
    p2 <- function(x) {formatC(x, format="f", digits=2)}
    p3 <- function(x) {formatC(x, format="f", digits=3)}
    p4 <- function(x) {formatC(x, format="f", digits=4)}
    p5 <- function(x) {formatC(x, format="f", digits=5)}
    
    logit <- function(p) log(1/(1/p-1))
    expit <- function(x) 1/(1/exp(x) + 1)
    inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
    is.even <- function(x){ x %% 2 == 0 } # function to identify odd maybe useful
    
    options(width=200)
    options(scipen=999)
    w=4  # line type
    ww=3 # line thickness
    wz=1 
    
    
    RR=.37  # used to limit correlations between variables
    sd1= 3  # for X covariates sd
    
    # links to Rdata objects uploaded to Git, these are pre run simulations see the save function
    # see line 1224 for the save function
    pp<- "https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/A%205000%20default%20settings%20theta%20log1.5%20-1.00%20-0.67%20-0.43.Rdata" # 5000 default log1.5 -1 -.67 -.43
    pp2<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/B%205000%20default%20settings%20theta%20log0.5%20-1.68%20-1.39%20%200.71.Rdata"
    pp3<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/C%205000%20default%20settings%20theta%20log2%20-3.46%20-1.05%20%201.15%20p1=.75.Rdata"
    pp4<-"https://github.com/eamonn2014/RCT-covariate-adjust-binary-response/raw/master/cov-adj-binary-response/D_10Ksims_5covariates_p1_0.12_theta_log1.3_covariates_-1.02%20_0.42_0.43_0.61%20_1.01_3_prog.Rdata"
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # load(url(pp))
    
    K1<-K<-3 
    Kp<-2
    p1<-0.35 
    theta<-log(1.5)
    pow <- .90 
    alpha<-5/100
    Fact<-3 
    simuls<-50 
    covar<-2 

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #power the examples
  post.odds <- p1/(1-p1) * exp(theta)
  post.prob <- post.odds/(1+post.odds)
  Po <- bsamsize(p1, post.prob, fraction=.5, alpha=alpha, power=pow)  # getting the sampe size
  
  MM <- N <-ceiling(Po[1][[1]])*2  # total will always be even
  bigN <- MM  
  N1 <- MM/2
  N2 <- N1
  
  # not used, as I want se for log odds ratio 
  # se. <- sqrt(  p1*(1-p1) /(N1) +    post.prob*(1- post.prob)/(N2) ) 
  
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SIMULATION CODE STARTS HERE
# here is code to simulate scenarios prognostic covariates, covariates unrelated to y, mix of pro and unrelated to y covariates, 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

N1 <-  N # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theta1 <- theta
# making up some beta coefficients, fixed for all simulations as it is outside loop
tx <- abs(-theta1)
b1 <- round(sort(runif(K1, -tx*Fact, tx*Fact)), digits=2) 

simfun <- function(N=N1, K=K1, a=log(p1/(1-p1)), sigma=sigma1, theta=(theta1), b=b1) {
  
  # randomi <- runif(N)  
  # we can select this, does not seem to have a big impact
  if (covar==1) {  
    X <- array(runif(N*K , -1,1), c(N,K))     # initially covars were uniform dist
  } else {
    X <- array(rnorm(N*K, 0, sd1), c(N,K))   
  }
  
  z <- sample(c(0,1), N, replace=T)            # treatment indicators
  
  # NEW CODE
  lp<-NULL
  lp = a+ X %*% b + theta*z 
  y <- ifelse(runif(N)  < plogis(lp), 1, 0)    # one liner RANDOM!!!
  lp<-NULL
  lp = a + theta*z 
  y2 <- ifelse(runif(N)  < plogis(lp), 1, 0)   # one liner RANDOM!!!
  lp<-NULL
  lp <- a+ X[,1:Kp] %*% b[1:Kp] + theta*z      # just have 3 variables associated with response
  y3 <- ifelse(runif(N)  < plogis(lp), 1, 0)   # one liner RANDOM!!!
  
  data.frame(X=X, y=y, z=z, y2=y2, y3=y3)
  
}

#https://stackoverflow.com/questions/5251507/how-to-succinctly-write-a-formula-with-many-variables-from-a-data-frame

# function to run the analyses
statfun <- function(d) {
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  zz <- glm(y~.-y2-y3, data=d,   family = "binomial")     ## adjusting for prognostic X, y2 is not included by use of the '-'
  f <-  summary(zz)
  
  zz1 <- glm(y~z, data=d,   family = "binomial")          ## not adjusting for prognostic X, only trt. indicator included
  f1 <-  summary(zz1)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  zz2 <- glm(y2~.-y-y3,data=d,   family = "binomial")     ## adjusting for X which are not prognostic
  f2 <-  summary(zz2)
  
  zz3 <- glm(y2~z, data=d,   family = "binomial")         ## not adjusting for X which are not prognostic
  f3 <-  summary(zz3)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  zz4 <- glm(y3~.-y-y2,data=d,   family = "binomial")     ## adjusting some X  are prognostic
  f4 <-  summary(zz4)
  
  zz5 <- glm(y3~z, data=d,   family = "binomial")         ## not adjusting when some X  are prognostic
  f5 <-  summary(zz5)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
  
  cbind(
    
    coef(f)["z", "Estimate"],
    coef(f)["z", "Std. Error"],
    
    coef(f1)["z", "Estimate"],
    coef(f1)["z", "Std. Error"],
    
    coef(f2)["z", "Estimate"],
    coef(f2)["z", "Std. Error"],
    
    coef(f3)["z", "Estimate"],
    coef(f3)["z", "Std. Error"], #8
    
    coef(f4)["z", "Estimate"],
    coef(f4)["z", "Std. Error"],
    
    coef(f5)["z", "Estimate"],
    coef(f5)["z", "Std. Error"],
    
    # collect p values
    coef(f)["z", "Pr(>|z|)"]  < alpha,  #13    
    coef(f1)["z", "Pr(>|z|)"] < alpha,
    coef(f2)["z", "Pr(>|z|)"] < alpha,
    coef(f3)["z", "Pr(>|z|)"] < alpha,
    coef(f4)["z", "Pr(>|z|)"] < alpha,
    coef(f5)["z", "Pr(>|z|)"] < alpha  , #18
    ## mse
    #mse for logistic regression https://rdrr.io/cran/dvmisc/man/get_mse.html
    sum(zz$residuals^2) / zz$df.residual ,
    sum(zz1$residuals^2) / zz1$df.residual, 
    sum(zz2$residuals^2) / zz2$df.residual ,
    sum(zz3$residuals^2) / zz3$df.residual ,
    sum(zz4$residuals^2) / zz4$df.residual ,
    sum(zz5$residuals^2) / zz5$df.residual ,
    
    (quantile( sum(zz$residuals^2) / zz$df.residual , .025)), #25
    (quantile( sum(zz$residuals^2) / zz$df.residual , .975)), 
    (quantile( sum(zz1$residuals^2) / zz1$df.residual , .025)), #25
    (quantile( sum(zz1$residuals^2) / zz1$df.residual , .975)), 
    (quantile( sum(zz2$residuals^2) / zz2$df.residual , .025)), #25
    (quantile( sum(zz2$residuals^2) / zz2$df.residual , .975)), 
    (quantile( sum(zz3$residuals^2) / zz3$df.residual , .025)), #25
    (quantile( sum(zz3$residuals^2) / zz3$df.residual , .975)), 
    (quantile( sum(zz4$residuals^2) / zz4$df.residual , .025)), #25
    (quantile( sum(zz4$residuals^2) / zz4$df.residual , .975)), 
    (quantile( sum(zz5$residuals^2) / zz5$df.residual , .025)), #25
    (quantile( sum(zz5$residuals^2) / zz5$df.residual , .975)), 
    
    f$aic,  #37
    f1$aic,
    f2$aic,
    f3$aic,
    f4$aic,
    f5$aic ,   #42
    
    # https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/
    1-(-f$deviance/2)/(-f$null.deviance/2),
    1-(-f1$deviance/2)/(-f1$null.deviance/2),
    1-(-f2$deviance/2)/(-f2$null.deviance/2),
    1-(-f3$deviance/2)/(-f3$null.deviance/2),
    1-(-f4$deviance/2)/(-f4$null.deviance/2),
    1-(-f5$deviance/2)/(-f5$null.deviance/2)
  )
}

library(plyr)
res <- raply(simuls, statfun(simfun())) # run the model many times
# summarize
result <- apply(res,2,mean)
q1.result <- apply(res,2, quantile, probs=c(0.025), na.rm=TRUE)
q2.result <- apply(res,2, quantile, probs=c(0.975), na.rm=TRUE)
# collect




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulation code 2 do the same as the first simulation code, but this time correlated covariates are created
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulate models many times collect estimate and SE
  
  simfun2 <- function(N=N1, K=K1, a=log(p1/(1-p1)), theta=(theta1), b=b1) {
    # randomi <- runif(N)  
    x <- Matrix(runif(K*K,-RR,RR), K)   # create a correlation matrix randomly , wont allow very high correlations
    
    A <- forceSymmetric(x)
    
    diag(A) <- 1
    
    M <- A
    
    M <- nearPD(M, conv.tol = 1e-7)$mat # default
    # Cholesky decomposition
    L = chol(M)
    nvars = dim(L)[1]
    
    # Random variables that follow an M correlation matrix
    r = t(L) %*% matrix(rnorm(nvars*N, 0,sd1), nrow=nvars, ncol=N)  #2,2
    r = t(r)
    
    r <- as.matrix(r)#
    rdata <- as.data.frame(r)
    XX<- as.matrix(rdata)
    z <- sample(c(0,1), N, replace=T)                # treatment indicator
    lp<-NULL
    lp <- a+ XX %*% b + theta*z     # betas created earlier
    y <- ifelse(runif(N)  < plogis(lp), 1, 0)   # one line
    
    data.frame(X=XX, y=y, z=z)
    
  }
  
  #https://stackoverflow.com/questions/5251507/how-to-succinctly-write-a-formula-with-many-variables-from-a-data-frame
  
  statfun2 <- function(d) {
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    zz <- glm(y~. , data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
    f <-  summary(zz)
    
    zz1 <- glm(y~z , data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
    f1 <-  summary(zz1)
    
    cbind(
      coef(f)["z", "Estimate"],
      coef(f)["z", "Std. Error"],
      coef(f1)["z", "Estimate"],
      coef(f1)["z", "Std. Error"], 
      
      coef(f)["z", "Pr(>|z|)"]  < alpha,   
      coef(f1)["z", "Pr(>|z|)"]  < alpha,
      
      ## mse
      sum(zz$residuals^2) / zz$df.residual ,
      sum(zz1$residuals^2) / zz1$df.residual, 
      
      (quantile( sum(zz$residuals^2) / zz$df.residual , .025)), #25
      (quantile( sum(zz$residuals^2) / zz$df.residual , .975)), 
      
      (quantile( sum(zz1$residuals^2) / zz1$df.residual , .025)), #25
      (quantile( sum(zz1$residuals^2) / zz1$df.residual , .975)), 
      
      f$aic,  #37
      f1$aic,
      
      # https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/
      1-(-f$deviance/2)/(-f$null.deviance/2),
      1-(-f1$deviance/2)/(-f1$null.deviance/2)
      
    )
    
  }
  
  library(plyr)
  res2 <- raply(simuls, statfun2(simfun2())) # run the model many times
  result2 <- apply(res2,2,mean)
  q1.result2 <- apply(res2,2, quantile, probs=c(0.025), na.rm=TRUE)
  q2.result2 <- apply(res2,2, quantile, probs=c(0.975), na.rm=TRUE)
  
 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# do the same as the first simulation code, but this time imbalanced covariates are created
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # simulate models many times collect estimate and SE
  
  simfun3<- function(N=N1, K=K1, a=log(p1/(1-p1)),  theta=(theta1), b=b1) {
    
    MM = N
    N2=MM/2
    N1=N2 
    
    if (covar==1) {  
      X1 <- array(runif(N1*K , -1,1), c(N1,K))  
      X2 <- array(runif(N2*K , -.8,1.2), c(N2,K))   ##imbalance compared to above
      XY <- X <- rbind(X1,X2)
    } else {
      X1 <- array(rnorm(N1*K, 0,  1), c(N1,K))  
      X2 <- array(rnorm(N2*K, .3, 1), c(N2,K))   ##imbalance compared to above
      XY <- X <- rbind(X1,X2)
    }
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #randomi <- runif(N) 
    z <- rep(0:1, c(N1,N2))  #assign 1 so we maintain shift in arms
    lp<-NULL
    lp = a+ X %*% b + theta*z 
    y <- ifelse(runif(N)  < plogis(lp), 1, 0)   # one liner RANDOM!!!
    lp<-NULL
    lp = a + theta*z 
    y2 <- ifelse(runif(N)  < plogis(lp), 1, 0)   # one liner RANDOM!!!
    
    data.frame(X=XY, y=y, z=z, y2=y2)
    
  }
  
  #https://stackoverflow.com/questions/5251507/how-to-succinctly-write-a-formula-with-many-variables-from-a-data-frame
  
  statfun3 <- function(d) {
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    zz <- glm(y~.-y2,  data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
    f <-  summary(zz)
    
    zz1 <- glm(y~z,  data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
    f1 <-  summary(zz1)
    
    zz2 <- glm(y2~.-y,  data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
    f2 <-  summary(zz2)
    
    zz3 <- glm(y2~z,  data=d,   family = "binomial")      ## adjusting for prognostic X, y2 is not included by use of the '-'
    f3 <-  summary(zz3)
    
    cbind(
      coef(f)["z", "Estimate"],
      coef(f)["z", "Std. Error"],
      
      coef(f1)["z", "Estimate"],
      coef(f1)["z", "Std. Error"],
      
      coef(f2)["z", "Estimate"],
      coef(f2)["z", "Std. Error"],
      
      coef(f3)["z", "Estimate"],
      coef(f3)["z", "Std. Error"], #8
      
      # collect p values for power
      coef(f)["z", "Pr(>|z|)"]  < alpha,  #9
      coef(f1)["z", "Pr(>|z|)"] < alpha,
      coef(f2)["z", "Pr(>|z|)"] < alpha,
      coef(f3)["z", "Pr(>|z|)"] < alpha,
      
      ## mse
      #mse for logistic regression https://rdrr.io/cran/dvmisc/man/get_mse.html
      sum(zz$residuals^2) / zz$df.residual ,
      sum(zz1$residuals^2) / zz1$df.residual, 
      sum(zz2$residuals^2) / zz2$df.residual ,
      sum(zz3$residuals^2) / zz3$df.residual ,
      
      
      (quantile( sum(zz$residuals^2) / zz$df.residual , .025)), #25
      (quantile( sum(zz$residuals^2) / zz$df.residual , .975)), 
      
      (quantile( sum(zz1$residuals^2) / zz1$df.residual , .025)), #25
      (quantile( sum(zz1$residuals^2) / zz1$df.residual , .975)), 
      
      (quantile( sum(zz2$residuals^2) / zz2$df.residual , .025)), #25
      (quantile( sum(zz2$residuals^2) / zz2$df.residual , .975)), 
      
      (quantile( sum(zz3$residuals^2) / zz3$df.residual , .025)), #25
      (quantile( sum(zz3$residuals^2) / zz3$df.residual , .975)), 
      
      f$aic,  #37
      f1$aic,
      f2$aic,
      f3$aic,
      
      # https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/
      1-(-f$deviance/2)/(-f$null.deviance/2),
      1-(-f1$deviance/2)/(-f1$null.deviance/2),
      1-(-f2$deviance/2)/(-f2$null.deviance/2),
      1-(-f3$deviance/2)/(-f3$null.deviance/2) 
      
    )
    
  }
  
  library(plyr)
  res3 <- raply(simuls, statfun3(simfun3())) # run the model many times
  # summarize
  result3 <- apply(res3,2,mean)
  q1.result3 <- apply(res3,2, quantile, probs=c(0.025), na.rm=TRUE)
  q2.result3 <- apply(res3,2, quantile, probs=c(0.975), na.rm=TRUE)
  # collect
 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulation plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# collect simulation trt effect estimates from simulation and plot!

  
  
  d1 <-  density(res[,1] )
  d2 <-  density(res[,3] )
  d3 <-  density(res[,5] )
  d4 <-  density(res[,7] )
  d5 <-  density(res[,9] )
  d6 <-  density(res[,11] )
  d7 <-  density(res2[,1] )
  d8 <-  density(res2[,3] )
  d9 <-   density(res3[,1] )
  d10 <-  density(res3[,3] )
  d11 <-  density(res3[,5] )
  d12 <-  density(res3[,7] )
  
  dz <- max(c(d1$y, d2$y, d3$y, d4$y, d5$y, d6$y, d7$y, d8$y  , d9$y, d10$y, d11$y, d12$y  ))
  dx <- range(c(d1$x,d2$x,  d3$x, d4$x, d5$x, d6$x, d7$x, d8$x   , d9$x, d10$x, d11$x, d12$x  ))
  
  
    plot((d1), xlim = dx, main=paste0("Density of treatment estimates, truth= ",p3(theta1),""), ylim=c(0,dz),lty=wz, lwd=ww,
         xlab="Treatment effect log odds",  
         ylab="Density")                           
    lines( (d2), col = "black", lty=w, lwd=ww)  
    lines( (d3), col = "red", lty=wz, lwd=ww)    
    lines( (d4), col = "red", lty=w, lwd=ww)          
    lines( (d5), col = "blue", lty=wz, lwd=ww)       
    lines( (d6), col = "blue", lty=w, lwd=ww)       
    lines( (d7), col = "purple", lty=wz, lwd=ww)       
    lines( (d8), col = "purple", lty=w, lwd=ww)       
    
    lines( (d9), col = "green", lty=wz, lwd=ww)       
    lines( (d10), col = "green", lty=w, lwd=ww)       
    lines( (d11), col = "grey", lty=wz, lwd=ww)       
    lines( (d12), col = "grey", lty=w, lwd=ww)  
    
 
  
  abline(v = theta1, col = "darkgrey")                
  legend("topright",       # Add legend to density
         legend = c(" adj. for true prognostic covariates", 
                    " not adj. for true prognostic covariates" ,
                    " adj. for covariates unrelated to outcome", 
                    " not adj. for covariates unrelated to outcome",
                    " adj. for mix of prognostic and unrelated to outcome", 
                    " not adj. mix of prognostic and unrelated to outcome", 
                    " adj. for correlated prognostic covariates", 
                    " not adj. for correlated prognostic covariates",
                    " adj. for imbalanced prognostic covariates", 
                    " not adj. for imbalanced prognostic covariates", 
                    " adj. for imbalanced covariates unrelated to outcome", 
                    " not adj. imbalanced covariates unrelated to outcome"
                    
         ),
         col = c("black", "black","red","red","blue", "blue", "purple", "purple", "green", "green", "grey", "grey"),
         lty = c(wz, w,wz,w,wz,w,wz,w,wz,w,wz,w)  ,lwd=ww
         , bty = "n", cex=1)
 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    

# Need to be able to calculate logistic regression se analytically base on p1 and OR?
# so for now we do simulation, not ideal

 
# 
#   N<- mcmc()$N
# 
#   ## function to simulate simple dataset #####################################################################
#   
#   fun.d<-function(nsample, drug.allocation, 
#                   alpha,  beta.drug,
#                   seed=NULL){ 
#     
#     if (!is.null(seed)) set.seed(seed)
#     
#     drug<- (rbinom(nsample, 1, prob =drug.allocation ))   
#     
#     Xmat <- model.matrix(~ drug )
#     beta.vec <- c(alpha,  beta.drug )
#     
#     lin.pred <- Xmat[,] %*% beta.vec                 # Value of lin.predictor
#     exp.p <- exp(lin.pred) / (1 + exp(lin.pred))     # Expected proportion
#     y <- rbinom(n = nsample, size = 1, prob = exp.p) # Add binomial noise
#     #y<- runif(nsample) <  exp.p                     # alternatively ads noise in this way
#     
#     d<-as.data.frame(cbind(y, drug))         # create a dataset
#     
#     return(d)
#     
#   }
#   
#  # function to pull out se######################################################################################
#   simfunc <- function(d) {
#     fit1 <- glm(y  ~ drug , d, family = binomial) 
#     c( summary(fit1)$coef["drug","Std. Error"] )   # collect se
#   }
#   
#  ################################################################################################################
#  out <- replicate(1000, simfunc(fun.d( nsample=N, drug.allocation=0.5,    
#                                       alpha=log(p1/(1-p1)),  # turn prob into odds and enter this here
#                                       beta.drug=(theta1))))
# 
# ################################################################################################################
#  se. <- mean(out)
# 
#    
# return(list(   
#   
#   se.=se. 
# ))

###############################################################################################################
# True standard error 
# not relying on simulation
# logistic regression se by hand for log odds ratio of drug effect based on p1 and OR from which N is calculated
############################################################################################################### 
 
  
 
  
  
  (post.odds <- p1/(1-p1) *  (theta1))          # calculate post odds
  (p2 <-post.prob <- post.odds/(post.odds+1))  # calculate post prob
  
  # calcualte N in each cell of 2by2 table
  cellA <-  (1-p1)*N/2
  cellC	<- p1*N/2
  cellB	<- (1-p2)*N/2
  cellD	<- p2*N/2
  
  # now se of log odds ratio
  se. <- log.odds.se. <- sqrt((1/cellA)+(1/cellB)+(1/cellC)+(1/cellD))
  
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SIMUALTION PLOT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# collect simulatio standard error estimates from simulation and plot!
 
 
  #######################
  # here we save simulation results
  # rename if important to keep
   #######################
  
  d1 <-  density(res[,2] )
  d2 <-  density(res[,4] )
  d3 <-  density(res[,6] )
  d4 <-  density(res[,8] )
  d5 <-  density(res[,10] )
  d6 <-  density(res[,12] )
  d7 <-  density(res2[,2] )
  d8 <-  density(res2[,4] )
  
  d9 <-   density(res3[,2] )
  d10 <-  density(res3[,4] )
  d11 <-  density(res3[,6] )
  d12 <-  density(res3[,8] )
  
  
  dz <- max(c(d1$y, d2$y, d3$y, d4$y, d5$y, d6$y, d7$y, d8$y  , d9$y, d10$y, d11$y, d12$y    ))
  dx <- range(c(d1$x,d2$x,  d3$x, d4$x, d5$x, d6$x, d7$x, d8$x   , d9$x, d10$x, d11$x, d12$x    ))
  
     
    plot( (d1), xlim = c(dx), main=paste0("Density of treatment standard error estimates, truth= ",p4(se.),""), ylim=c(0,dz),lty=wz, lwd=ww,
          xlab="Standard error of log odds trt effect",  
          ylab="Density")  
    lines( (d2), col = "black", lty=w, lwd=ww)  
    lines( (d3), col = "red", lty=wz, lwd=ww)    
    lines( (d4), col = "red", lty=w, lwd=ww)          
    lines( (d5), col = "blue", lty=wz, lwd=ww)       
    lines( (d6), col = "blue", lty=w, lwd=ww)       
    lines( (d7), col = "purple", lty=wz, lwd=ww)       
    lines( (d8), col = "purple", lty=w, lwd=ww)       
    
    lines( (d9), col = "green", lty=wz, lwd=ww)       
    lines( (d10), col = "green", lty=w, lwd=ww)       
    lines( (d11), col = "grey", lty=wz, lwd=ww)       
    lines( (d12), col = "grey", lty=w, lwd=ww)  
    
  
  abline(v = se., col = "darkgrey")   
  legend("topright",           # Add legend to density
         legend = c(" adj. for true prognostic covariates", 
                    " not adj. for true prognostic covariates" ,
                    " adj. for covariates unrelated to outcome", 
                    " not adj. for covariates unrelated to outcome",
                    " adj. for mix of prognostic and unrelated to outcome", 
                    " not adj. mix of prognostic and unrelated to outcome", 
                    " adj. for correlated prognostic covariates", 
                    " not adj. for correlated prognostic covariates",
                    " adj. for imbalanced prognostic covariates", 
                    " not adj. for imbalanced prognostic covariates", 
                    " adj. for imbalanced covariates unrelated to outcome", 
                    " not adj. imbalanced covariates unrelated to outcome"
                    
         ),
         col = c("black", "black","red","red","blue", "blue", "purple", "purple", "green", "green", "grey", "grey"),
         lty = c(wz, w,wz,w,wz,w,wz,w,wz,w,wz,w) ,lwd=ww
         , bty = "n", cex=1)
 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  p2 <- function(x) {formatC(x, format="f", digits=2)}  
  zz <- rbind(
    (c( p4(result[1])   ,     p2(q1.result[1])  ,  p2(q2.result[1])   , p4(result[2] ) ,  p4(result[13] ) ,  p4(result[1] -theta1) ,      p4(result[37] )    ,  p4(result[43] )         )) ,
    (c( p4(result[3])   ,     p2(q1.result[3]) ,   p2(q2.result[3])   , p4(result[4] ) ,  p4(result[14] ) ,  p4(result[3] -theta1) ,      p4(result[38] )    ,  p4(result[44] )        )) ,
    (c( p4(result[5])   ,     p2(q1.result[5]) ,   p2(q2.result[5])   , p4(result[6] ) ,  p4(result[15] ) ,  p4(result[5] -theta1) ,      p4(result[39] )    ,  p4(result[45] )        )) ,
    (c( p4(result[7])   ,     p2(q1.result[7]) ,   p2(q2.result[7])   , p4(result[8] ) ,  p4(result[16] ) ,  p4(result[7] -theta1) ,      p4(result[40] )    ,  p4(result[46] )         )) ,
    (c( p4(result[9])   ,     p2(q1.result[9]) ,   p2(q2.result[9])   , p4(result[10] ) , p4(result[17] ) ,  p4(result[9] -theta1) ,      p4(result[41] )    ,  p4(result[47] )           )) ,
    (c( p4(result[11])  ,     p2(q1.result[11]) ,  p2(q2.result[11])  , p4(result[12] ) , p4(result[18] ) ,  p4(result[11] -theta1) ,      p4(result[42] )    ,  p4(result[48] )         )) ,
    (c( p4(result2[1])  ,     p2(q1.result2[1]),   p2(q2.result2[1])  , p4(result2[2] ) , p4(result2[5] ) ,  p4(result2[1] -theta1) ,      p4(result2[13] )   ,  p4(result2[15] )         )) ,
    (c( p4(result2[3])  ,     p2(q1.result2[3])  , p2(q2.result2[3])  , p4(result2[4] ) , p4(result2[6] ) ,  p4(result2[3] -theta1) ,      p4(result2[14] )   ,  p4(result2[16] )          )),
    (c( p4(result3[1])  ,     p2(q1.result3[1])  , p2(q2.result3[1])   , p4(result3[2] ) ,  p4(result3[9] ) ,  p4(result3[1] -theta1) , p4(result3[25] )   ,  p4(result3[29] )         )) ,
    (c( p4(result3[3])  ,     p2(q1.result3[3]) ,  p2(q2.result3[3])   , p4(result3[4] ) ,  p4(result3[10] ) ,  p4(result3[3] -theta1) ,p4(result3[26] )   ,  p4(result3[30] )         )) ,
    (c( p4(result3[5])  ,     p2(q1.result3[5]) ,  p2(q2.result3[5])   , p4(result3[6] ) ,  p4(result3[11] ) ,  p4(result3[5] -theta1) , p4(result3[27] )  ,  p4(result3[31] )         )) ,
    (c( p4(result3[7])  ,     p2(q1.result3[7]) ,  p2(q2.result3[7])   , p4(result3[8] ) ,  p4(result3[12] ) ,  p4(result3[7] -theta1) , p4(result3[28] )  ,  p4(result3[32]  )         )) 
  ) 
  
  zz <- as.data.frame(zz)
  
  colnames(zz) <- c("Mean  ", "Lower 95%CI", "Upper 95%CI", "Stand.error", "Power ","B" , "sigma","R2")
  
  zz <- data.frame(lapply(zz, function(x) as.numeric(as.character(x))))
  zz <- as.data.frame(zz)
  rownames(zz)<- c(
    " adj. for true prognostic covariates", 
    " not adj. for true prognostic covariates" ,
    " adj. for non prognostic covariates", 
    " not adj. for non prognostic covariates",
    " adj. for some non prognostic covariates", 
    " not adj. when some prognostic covariates", 
    " adj. for correlated prognostic covariates", 
    " not adj. for correlated prognostic covariates",
    " adj. for imbalanced prognostic covariates",
    " not adj. for imbalanced prognostic covariates",  
    " adj. for non prognostic imbalanced covariates",
    " not adj. for non prognostic imbalanced covariates"
  )
  zz <- zz[order(zz$B),]
  
  
  colnames(zz) <- c("Mean  ", "Lower 95%CI", "Upper 95%CI", "Std.error", "Power ","bias" , "AIC","McFadden's R2")
  zz <- zz[order(zz$AIC),]
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # simulation 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
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
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # functions to fit model and pull out a stat
  # lrtest
  simfunc <- function(d) {
    fit1 <- glm(y  ~ drug , d, family = binomial) 
    fit2 <- glm(y  ~ 1, d, family = binomial) 
    c(anova(fit1,fit2, test='Chisq')[2,5] )
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # wald test
  simfunc <- function(d) {
    fit1 <- glm(y  ~ drug , d, family = binomial) 
    c( summary(fit1)$coef["drug","Pr(>|z|)"] )
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  NN <- 300
  
  pp1 <- .25
  
  pp2 <- .35   # odds ratio 
  
  allocation <- .5

  orx <- "yes"  # what do we want pp2 to be odds ratio or baseline proportion
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (orx %in% "yes") {
    pp3 <- pp1 * pp2/(1 - pp1 + pp1 * pp2)
    or <- pp2
  } else  {
    pp3 <- pp2
    or <- pp2/(1-pp2) / (   pp1/(1-pp1) )
  }  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pp2 <- pp3
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    out <- replicate(499, simfunc(fun.d( nsample=NN,     
                                         drug.allocation=allocation,  
                                         alpha=log(odds),  
                                         beta.drug=log(or))))
    
    pow <-  mean( out < 0.05, na.rm=TRUE )   
    pow
    
    # frank harrell
    Hmisc::bpower(p1=pp1, odds.ratio=c(or ),  alpha=c(0.05) , n1=NN*(1-allocation), n2=NN*(allocation) ) 
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NN <- 300
    
    pp1 <- .25
    
    pp2 <- .35   # odds ratio 
 
    allocation <- .5
    
    orx <- "no"  # what do we want pp2 to be odds ratio or baseline proportion
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if (orx %in% "yes") {
      pp3 <- pp1 * pp2/(1 - pp1 + pp1 * pp2)
      or <- pp2
    } else  {
      pp3 <- pp2
      or <- pp2/(1-pp2) / (   pp1/(1-pp1) )
    }  
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pp2 <- pp3
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    out <- replicate(499, simfunc(fun.d( nsample=NN,     
                                         drug.allocation=allocation,  
                                         alpha=log(odds),  
                                         beta.drug=log(or))))
    
    pow <-  mean( out < 0.05, na.rm=TRUE )   
    pow
    
    # frank harrell
    Hmisc::bpower(p1=pp1, odds.ratio=c(or ),  alpha=c(0.05) , n1=NN*(1-allocation), n2=NN*(allocation) ) 
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  