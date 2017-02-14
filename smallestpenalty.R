warsmallest_penalty <- function(w,lambda,mu,t,maxeval=5000,...){
  p <- prob_s_tau(t,lambda,mu)
  smallest_penalty_obj<- function(x){-min(x)}
#  browser()
  EV <- function(h,w_prime,tau)
  {
    delta_A^(t-1)*sum(p[,tau+1]*(h + delta_A *w_prime)) - ifelse(tau==t,0,a*sum(delta_A^(tau:(t-1))))
  }
  
  heq <- function(x)
  {
    h <- x[1:(t+1)]
    w_prime <- x[(t+2):(2*t+2)]
    
    heq<- rep(NA,1)
    heq[1] <-  (w - EV(h,w_prime,tau=0))
#    if(t>1) for(i in 1:(t-1)) {
#      heq[1+i] <- h[i] - h[i+1]
#      heq[t+i] <- w_prime[i] - w_prime[i+1]
#    }
    heq
  }

  hin <- function(x) ## Inequality constraints ##
  {
    h <- x[1:(t+1)]
    w_prime <- x[(t+2):(2*t+2)]
    
    hin <- rep(NA,2*t+1)
    hin[1:(t+1)] <- w_prime
    for(tau in 1:t) hin[(t+1)+tau] <- (EV(h,w_prime,0)-EV(h,w_prime,tau)) ## ICC's for each length on deviation #
    
    hin
  }
  
  effort_high <- nloptr::auglag( x0=rep(0,2*(t+1)), 
                                 fn= smallest_penalty_obj, 
                                 hin = hin,
                                 heq = heq,
                                 nl.info = FALSE,
                                 control=list(xtol_rel = 1e-8, maxeval = maxeval))
    effort_high
}


lambda_seq <- seq(mu+0.1,1,0.1)
penalty <- matrix(ncol=3,nrow=length(lambda_seq))
for(i in 1:length(lambda_seq)){
  for(j in 1:3){
    penalty[i,j] <- smallest_penalty(0,lambda_seq[i],mu,j,maxeval=5000)$value
  }
}
plot(lambda_seq,penalty[,1],type="b")
lines(lambda_seq,penalty[,2]/2)
lines(lambda_seq,penalty[,3]/3)