

w_prime_s0 <- Schumaker(W,sol_t2_lambda0.9$Policy[,4])$Spline
w_prime_s1 <- Schumaker(W,sol_t2_lambda0.9$Policy[,5])$Spline
w_prime_s2 <- Schumaker(W,sol_t2_lambda0.9$Policy[,6])$Spline

### 45 degree line plot ####
par(mfrow=c(1,1))
plot(W[-(13:15)],W[-(13:15)],type="l")
lines(W[-(13:15)],sol_t1_lambda0.7$Policy[-(13:15),3])
lines(W[-(13:15)],sol_t1_lambda0.7$Policy[-(13:15),4])
lines(W,sol_t2_lambda0.5$Policy[,6])


###

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