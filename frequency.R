
library("ggplot2")
library("reshape2")

transition_plot <- function(lambda,t,W.=W,max_w=Inf) {
  
  W. <- W.[which(W<max_w)]
  
  ### 45 degree line plot ####
  
  sol <- eval(as.name(paste("sol_t",t,"_lambda",lambda,sep="")))
  
  df <- as.data.frame(cbind(W.,W.,sol$Policy[,(t+2):(2*t + 2)]))
  
  if(t==1) {colnames(df)<- c("W","45 degree","bad","good")
  } else {colnames(df) <- c("W","45 degree","bad","medium","good")}
  
  df <- melt(df ,  id.vars = 'W', variable.name = 'series')
  
  plot<- ggplot(df,aes(W, value))+geom_line(aes(colour=series))+ggtitle("Transition functions")
  
  plot
}

###

simulation <- function(lambda,t,W,iter=1000,no.seeds=10){

  sol <- eval(as.name(paste("sol_t",t,"_lambda",lambda,sep="")))
  sol$Policy<-round(sol$Policy,3)
  if(t==1) { 
  w_prime_s0 <- approxfun(W,sol$Policy[,3])
  w_prime_s1 <- approxfun(W,sol$Policy[,4])
  } else{
    w_prime_s0 <- approxfun(W,sol$Policy[,4])
    w_prime_s1 <- approxfun(W,sol$Policy[,5])
    w_prime_s2 <- approxfun(W,sol$Policy[,6])
  } 
  w_prime <- list(s0 = w_prime_s0,
                s1 = w_prime_s1,
                s2 = (ifelse(t>1,w_prime_s2,NA)))

  x <-  matrix(rep(0,iter*no.seeds),ncol=no.seeds);
  x[1,] <- runif(no.seeds,min=0,max=U(pi))
    for(i in 2:nrow(x)){
      outcome <- sum(as.integer(runif(t)<lambda))
      x[i,] <- w_prime[[outcome+1]](x[i-1,])
    }
  x
}
  
#dens <- logspline(x)
dens <- density(x)
plot(dens)
hist(x)
densit<- approxfun(dens$x,dens$y)
integrate(function(x)(densit(x)*sol$f(x)),min(dens$x),max(dens$x))
