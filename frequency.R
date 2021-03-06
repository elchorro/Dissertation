

transition_plot <- function(lambda,t,W.=W,max_w=Inf) {
  W. <- W.[which(W<max_w)]
  
  ### 45 degree line plot ####
  
  sol <- eval(as.name(paste("sol_t",t,"_lambda",lambda,sep="")))
  
  Policy <- sol$Policy[which(W<max_w),(t+2):(2*t + 2)]
  
  df <- as.data.frame(cbind(W.,W.,Policy))
  
  if(t==1) {colnames(df)<- c("W","45 degree","bad","good")
  } else {colnames(df) <- c("W","45 degree","bad","medium","good")}
  
  df <- melt(df ,  id.vars = 'W', variable.name = 'Signal')
  
  plot.subtitle <- paste("Lambda= ",lambda,"; t= ",t,".",sep="")
  plot<- ggplot(df,aes(W, value))+geom_line(aes(colour=Signal))+ggtitle(bquote(atop(.("Transition functions"), atop(italic(.(plot.subtitle)), ""))))
  
  plot
}
consumption_plot <- function(lambda,t,W.=W,max_w=Inf) {
  W. <- W.[which(W<max_w)]
  
  ### 45 degree line plot ####
  
  sol <- eval(as.name(paste("sol_t",t,"_lambda",lambda,sep="")))
  
  consumption <- cost(sol$Policy[which(W<max_w),1:(t + 1)])
  df <- as.data.frame(cbind(W.,consumption))
  
  if(t==1) {colnames(df)<- c("W","bad","good")
  } else {colnames(df) <- c("W","bad","medium","good")}
  
  df <- melt(df ,  id.vars = 'W', variable.name = 'Signal')
  
  plot.subtitle <- paste("Lambda= ",lambda,"; t= ",t,".",sep="")
  plot<- ggplot(df,aes(W, value))+geom_line(aes(colour=Signal))+ggtitle(bquote(atop(.("Consumption functions"), atop(italic(.(plot.subtitle)), ""))))
  
  plot
}
utility_plot <- function(lambda,t,W.=W,max_w=Inf) {
  W. <- W.[which(W<max_w)]
  
  ### 45 degree line plot ####
  
  sol <- eval(as.name(paste("sol_t",t,"_lambda",lambda,sep="")))
  
  consumption <- sol$Policy[which(W<max_w),1:(t + 1)]-a
  df <- as.data.frame(cbind(W.,consumption))
  
  if(t==1) {colnames(df)<- c("W","bad","good")
  } else {colnames(df) <- c("W","bad","medium","good")}
  
  df <- melt(df ,  id.vars = 'W', variable.name = 'Signal')
  
  plot.subtitle <- paste("Lambda= ",lambda,"; t= ",t,";a =",a,sep="")
  plot<- ggplot(df,aes(W, value))+geom_line(aes(colour=Signal))+ggtitle(bquote(atop(.("Utility provided (net effort cost)"), atop(italic(.(plot.subtitle)), ""))))
  
  plot
}
costf_plot <- function(lambda,t,W.=W,max_w=Inf) {
  W. <- W.[which(W<max_w)]
  
  ### 45 degree line plot ####
  
  sol <- eval(as.name(paste("sol_t",t,"_lambda",lambda,sep="")))
  
  values <- sol$values[which(W<max_w)]
  
  df <- as.data.frame(cbind(W.,values))
  
  colnames(df)<- c("W","Expected_cost ")
  
  df <- melt(df ,  id.vars = 'W', variable.name = '')
  
  plot.subtitle <- paste("Lambda= ",lambda,"; t= ",t,".",sep="")
  plot<- ggplot(df,aes(W, value))+geom_line(aes(y=value))+ggtitle(bquote(atop(.("Expected cost"), atop(italic(.(plot.subtitle)), ""))))
  
  plot
}

simulation <- function(lambda,t,W,iter=1000,no.seeds=10){

  sol <- eval(as.name(paste("sol_t",t,"_lambda",lambda,sep="")))
  
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
#dens <- density(x)
#plot(dens)
#hist(x)
#densit<- approxfun(dens$x,dens$y)
#integrate(function(x)(densit(x)*sol$f(x)),min(dens$x),max(dens$x))
