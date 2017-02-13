library("nloptr", lib.loc="~/R/win-library/3.3")
library("spsi")
library("schumaker")
library("alabama")
library("tictoc")
library("ConSpline")

a <- 1 # cost of effort 

pi <- 5 # output to principal

delta_A <- 0.8; delta_P <- 0.8 # discount rates

y <- function(e){ ifelse(e==1,pi,0) }

U <- function(c) c^(1/2)
cost <- function(u) u^2 #cost of providing giben utility, i.e. inverse of U 

bmax <- 100*pi
wmax <- U(30*pi)



## state space

## probability that signal i occurrs given tau deviations is p[i+1,tau+1], since "0" index not available
### needed for ICC ###

prob_s_tau <- function(t,lambda=lambda,mu=mu){ 
  p <- matrix(ncol=t+1,nrow=t+1)
  for(tau in 0:t){
    for(i in 1:(t+1))
      p[i,tau+1] <- sum(dbinom(0:tau,tau,mu)*dbinom(i-1-(0:tau),t-tau,lambda))
  }
  p
}


##### Bellman Operator #######

## Objective and constraints in Bellman operator

# x is a vector of length 2t+2. x[1:t+1] are current utilites, while x[t+2:2t+2] are promised 
# utilities under signals 1 and 2 respectively.


Bellman_P_alabama_multi <- function(w, Value_P, t=1, max_penal = 5*a, initial=NA,maxeval=5000,parameters=NULL){
  

  if(is.list(parameters)) for (i in 1:length(parameters)) assign(names(parameters)[i],parameters[[i]])  
  
  if(is.na(initial)) {
    initial[1:(t+1)] <- ((1:(t+1))^4)*(U(bmax)-0.0001)/((t+1)^4) 
    initial[(t+2):(2*(t+1))] <- ((1:(t+1))^4)*(wmax-0.00001)/((t+1)^4)
  }
  p <- prob_s_tau(t,lambda,mu)

  EV <- function(h,w_prime,tau) ## expected utility net effort cost
  {
    delta_A^(t-1)*sum(p[,tau+1]*(h + delta_A *w_prime)) - ifelse(tau==t,0,a*sum(delta_A^(tau:(t-1))))
  }
  
  
  objective <- function(x) 
  {
    h <- x[1:(t+1)]
    w_prime <- x[(t+2):(2*t+2)]
    obj <- delta_P^(t-1)*sum(p[,1]*(cost(h) + delta_P * Value_P(w_prime))) ## note p[,1] is vector of probabilities with no deviation
    obj
               
  }
  
  heq <- function(x)
  {
    h <- x[1:(t+1)]
    w_prime <- x[(t+2):(2*t+2)]
    
    w - EV(h,w_prime,0)
  }
  
  hin <- function(x) ## Inequality constraints ##
  {
    h <- x[1:(t+1)]
    w_prime <- x[(t+2):(2*t+2)]

    hin <- rep(NA,5*t+4)
    hin[1:(4*t+4)] <- c(h+max_penal, U(bmax)-h, w_prime, wmax-w_prime)
    for(tau in 1:t) hin[(4*t+4)+tau] <- (EV(h,w_prime,0)-EV(h,w_prime,tau)) ## ICC's for each length on deviation ##
    hin
  }

#  effort_high <- alabama::constrOptim.nl(par=initial, 
#                                         fn= objective, 
#                                         hin = hin,
#                                         heq = heq,
#                                         control.outer = list(itmax=maxeval,trace=F)
#  )
 
  effort_high <- nloptr::auglag( x0=initial, 
                  fn= objective,
                  hin = hin,
                  heq = heq,
                  nl.info = FALSE,
                  control=list(xtol_rel = 1e-8, maxeval = maxeval)) 
  
  try(if(min(hin(effort_high$par))< -0.0001) 
    warning(paste("ICC doesn't hold:",min(hin(effort_high$par)),", w=",w)))
  
  try(if(abs(heq(effort_high$par))>0.0001) 
    warning(paste("Promise-keeping doesn't hold:",abs(heq(effort_high$par)), ", w=",w)))
  effort_high
}

#### Iterating value function

value_new_alabama_multi <- function(fun=function(x)(x^2),iter=10, space=W,t=1,...){
  Policy <- matrix(ncol=2*(t+1),nrow=length(space))
  
  for(k in 1:iter){
    
    values <- rep(NA,length(W))
    
    for(i in 1:length(space)){
      optimal <- Bellman_P_alabama_multi(W[i],fun,t,...)
      values[i] <- optimal$value
      Policy[i,] <- optimal$par
    }
#    Values <- conspline(values,W,3)$muhat
    fun <- Schumaker(W,values)$Spline
    print(paste("Iteration ",k,"/", iter ," completed",sep=""))
  }
  
  parameters <- list("a"=a,"delta_A"=delta_A,"delta_P"=delta_P,"lambda"=lambda,"mu=mu")
  
  return(list("f"=fun,"values"=values,"Policy"=Policy))
}

###################################################################################
####################          Analysis            #################################
###################################################################################

for(lambda in c(0.5,0.7)){
  
  W <- seq(0,wmax,length.out=30)
  par(mfrow=c(1,2))
  
  for(t in 1:2){
    if(exists(paste("sol_t",t,"_lamda",lambda,sep=""))) f<- eval(as.name(paste("sol_t",t,"_lamda",lambda,sep="")))$f
    else f <- function(x)(exp(x))
    
    
    sol <- value_new_alabama_multi(fun=f,iter=10,t=t,maxeval=5000)
    assign(paste("sol_t",t,"_lamda",lambda,sep = ""),c(sol,list(parameters=list(a=a,delta_A=delta_A,delta_P=delta_P,lambda=lambda,mu=mu,W=W))))
    plot(W,sol$values)
    title(main = paste("lambda=",lambda,"t=",t))
  }
  
}
