##Kim 2006
##power calculation for case-cohort when the event is NOT rare
ccpw_freq <- function(alpha=0.95,hr,p1,N,pfail,q){
  #hr assumpted hazard ratio
  theta <- log(hr)
  #p1 proportion of the population in group1
  #p2 proportion of the population in group2
  p2 <- 1-p1
  #N total number of subject in the full cohort
  #n total number of subject in the subcohort
  #q sampling fraction of the subcohort
  n <- N*q
  #pfail proportion of failures in the full cohort
  Zalpha <- qnorm(alpha,lower.tail = F)
  
  #lambda can be found by solving
  #1-pd=lambda^1*(1-exp(-lambda))
  pfail_function <- function(lambda,pd){
    1-lambda^-1*(1-exp(-lambda))-pd
  }
  lambda <- uniroot(pfail_function,interval=c(-1,5),pd=pfail)$root
  
  #define A
  A <- exp(-lambda)+2*pfail-1
  
  #calculate power
  pwr <- pnorm(
    Zalpha+sqrt(n)*theta*
      sqrt((p1*p2*pfail)/
             (q+(1-q)*2*A/pfail))
  )
  return(pwr)
}
