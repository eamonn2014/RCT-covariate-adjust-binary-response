
OR<-2
prior.p<-seq(0.1,0.9, by=0.1)
prior.odds<-prior.p/(1-prior.p)
post.odds<-prior.odds*OR
post.p<-(post.odds/(1+post.odds))

f<-cbind(prior.p, prior.odds, OR, post.odds, post.p)
colnames(f)<-c("prior prob.","prior odds","odds ratio","post odds","predicted prob.")
print(f, digits=2)


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