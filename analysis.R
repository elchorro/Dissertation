###################################################################################
####################          Analysis            #################################
###################################################################################

a <- 1 # cost of effort 
mu <- 0.2 # probability of success given lowe effort
pi <- 5 # output to principal

delta_A <- 0.9; delta_P <- 0.9 # discount rates


y <- function(e){ ifelse(e==1,pi,0) }

U <- function(c) log(c+1)

cost <- function(u){
   ifelse(u>=0, exp(u)-1,0)}
  
#cost of providing giben utility, i.e. inverse of U 

bmax <- 100000*pi
wmax <- U(8*pi)




transitionplots <- list()

W <- seq(0,wmax,length.out=10)

##

for(lambda in c(0.5,0.7)){

  par(mfrow=c(1,2))
  
  for(t in 1:2){
    if(exists(paste("sol_t",t,"_lambda",lambda,sep=""))) f<- eval(as.name(paste("sol_t",t,"_lambda",lambda,sep="")))$f
    else f <- function(x)(x^2)
    
    
    sol <- value_new_alabama_multi(fun=f,iter=4,t=t,maxeval=7000,localsolver.= "MMA")
    assign(paste("sol_t",t,"_lambda",lambda,sep = ""),c(sol,list(parameters=list(a=a,delta_A=delta_A,delta_P=delta_P,lambda=lambda,mu=mu,W=W))))
    plot(W,sol$values)
    title(main = paste("lambda=",lambda,"t=",t))
    
    
  }
  
  
}


