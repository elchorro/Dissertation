###################################################################################
####################          Analysis            #################################
###################################################################################

a <- 1 # cost of effort 
mu <- 0.2 # probability of success given lowe effort
pi <- 5 # output to principal

delta_A <- 0.7; delta_P <- 0.7 # discount rates

y <- function(e){ ifelse(e==1,pi,0) }

U <- function(c) log(c+1)

cost <- function(u){
   ifelse(u>=0, exp(u)-1,0)}
  
#cost of providing giben utility, i.e. inverse of U 

bmax <- 1000*pi
wmax <- U(30*pi)




transitionplots <- list()

for(lambda in c(0.5,0.7)){
  
  W <- seq(0,wmax,length.out=20)
  par(mfrow=c(1,2))
  
  for(t in 1:2){
    if(exists(paste("sol_t",t,"_lambda",lambda,sep=""))) f<- eval(as.name(paste("sol_t",t,"_lambda",lambda,sep="")))$f
    else f <- function(x)(exp(x))
    
    
    sol <- value_new_alabama_multi(fun=f,iter=5,t=t,maxeval=5000)
    assign(paste("sol_t",t,"_lambda",lambda,sep = ""),c(sol,list(parameters=list(a=a,delta_A=delta_A,delta_P=delta_P,lambda=lambda,mu=mu,W=W))))
    plot(W,sol$values)
    title(main = paste("lambda=",lambda,"t=",t))
    
    transitionplots<- list(transitionplots, list(lambda=lambda,t=t,plot= transition_plot(lambda,t,W,7)))
    
  }
  
  
}


