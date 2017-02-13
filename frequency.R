
library("ggplot2")
library("reshape2")

transition_plot <- function(lambda,t,W)
{
  l<- length(W)
  ### 45 degree line plot ####
  par(mfrow=c(1,1))
  
  
  colours <- c("#999999", "#E69F00", "#56B4E9")
  
  sol <- eval(as.name(paste("sol_t",t,"_lambda",lambda,sep="")))
  
  df <- as.data.frame(cbind(W,sol$Policy[,4:6]))
  df <- melt(df ,  id.vars = 'W', variable.name = 'series')
  plot<- ggplot(df,aes(W, series))+geom_line(aes(colour=series))
  
  }
  plot
}

###


w_prime_s0 <- Schumaker(W,sol_t2_lambda0.9$Policy[,4])$Spline
w_prime_s1 <- Schumaker(W,sol_t2_lambda0.9$Policy[,5])$Spline
w_prime_s2 <- Schumaker(W,sol_t2_lambda0.9$Policy[,6])$Spline


sol <- eval(as.name(paste("sol_t",t,"_lambda",lambda,sep="")))

w_prime <- list(s0 = Schumaker(W,sol$Policy[,t+1])$Spline,
                s1 = Schumaker(W,sol$Policy[,t+2])$Spline,
                s2 = (ifelse(t>1,Schumaker(W,sol$Policy[,t+3])$Spline,NA)))

no_periods <- 1000

x <-  rep(NA,no_periods)
x[1] <- wmin
for(i in 2:length(x)){
    outcome <- sum(as.integer(runif(t)<lambda))
    x[i] <- w_prime[[outcome+1]](x[i-1])
}

#dens <- logspline(x)
dens <- density(x)
plot(dens)
hist(x)
densit<- approxfun(dens$x,dens$y)
integrate(function(x)(densit(x)*sol$f(x)),min(dens$x),max(dens$x))