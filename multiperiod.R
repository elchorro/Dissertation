library("nloptr", lib.loc="~/R/win-library/3.3")
library("spsi")
library("schumaker")
library("alabama")
library("tictoc")
library("ConSpline")  ###vv

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


Bellman_P_alabama_multi <- function(w, Value_P, t=1, max_penal = 5*a, initial=NA,iter=5000,parameters=NULL,wmax.=wmax,hmax=U(bmax)){
  

  if(is.list(parameters)) for (i in 1:length(parameters)) assign(names(parameters)[i],parameters[[i]])  
  
  if(is.na(initial)) {
    initial[1:(t+1)] <- ((1:(t+1))^4)*(hmax-0.0001)/((t+1)^4) 
    initial[(t+2):(2*(t+1))] <- ((1:(t+1))^4)*(wmax.-0.00001)/((t+1)^4)
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
    hin[1:(4*t+4)] <- c(h+max_penal, hmax-h, w_prime, wmax.-w_prime)
    for(tau in 1:t) hin[(4*t+4)+tau] <- (EV(h,w_prime,0)-EV(h,w_prime,tau)) ## ICC's for each length on deviation ##
    hin
  }

#  effort_high <- alabama::constrOptim.nl(par=initial, 
#                                         fn= objective, 
#                                         hin = hin,
#                                         heq = heq,
#                                         control.outer = list(itmax=iter,trace=F)
#  )
 
  effort_high <- nloptr::auglag( x0=initial, 
                  fn= objective,
                  hin = hin,
                  heq = heq,
                  nl.info = FALSE,
                  control=list(xtol_rel = 1e-8, maxeval = iter)) 
  

  try(if(min(head(hin(effort_high$par),4*(t+1))) < -0.0001) 
    warning(paste("Utility bound is binding:",min(hin(effort_high$par)),",at w=",w)))
 
  try(if(min(tail(hin(effort_high$par),t+1)) < -0.0001) 
    warning(paste("ICC is binding:",min(hin(effort_high$par)),",at w=",w)))
   
  try(if(abs(heq(effort_high$par))>0.0001) 
    warning(paste("Promise-keeping doesn't hold:",abs(heq(effort_high$par)), ",at  w=",w)))
  effort_high
}

#### Iterating value function

value_new_alabama_multi <- function(fun=function(x)(x^2),iter=10, space=W,t=1,maxeval=5000,...){
  Policy <- matrix(ncol=2*(t+1),nrow=length(space))
  
  for(k in 1:iter){
    
    values <- rep(NA,length(W))
    
    for(i in 1:length(space)){
      if(i==1) {
        optimal <- Bellman_P_alabama_multi(W[i],fun,t,iter=min(maxeval,10000),...)
      } else {optimal <- Bellman_P_alabama_multi(W[i],fun,t,iter=maxeval,...) }
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

