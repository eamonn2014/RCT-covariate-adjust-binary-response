
OR<-5
prior.p <- c(0.01,seq(0.05,0.95, by=0.05), .99)
prior.odds <- prior.p/(1-prior.p)
post.odds <- prior.odds*OR
post.p <- (post.odds/(1+post.odds))

f<-cbind(prior.p, prior.odds, OR, post.odds, post.p)
colnames(f)<-c("prior prob.","prior odds","odds ratio","post odds","predicted prob.")
print(f, digits=2)

par(pty="s")
plot(f[,1], f[,5], xlab="prior prob.", ylab='post.prob.', type = "l" , ylim=c(0,1) , xlim=c(0,1))
abline(0,1, lty=2)#

#lines( f[,2], f[,4]       )
#if the OR is 1.5 for a 0.1 unit change then...
#for a x unit change 'want'

or<-1.75 ;unit<-0.1
want<-0.2
ratio<-want/unit

or<-or^ratio
or

 

post.odds <- exp(-3)/(1-exp(-3))*exp(1)

post.p<-(post.odds/(1+post.odds))

require(Hmisc)
#showing equivalence in p1 OR and p1 p2
bpower(p1=.1, odds.ratio=2, n=1000, alpha=c(.01,.05))
bpower(p1=.1, p2=0.1818182, n=1000, alpha=c(.01,.05))