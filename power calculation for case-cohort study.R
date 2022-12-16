##Cai and Zeng (2004)
##power calculation for case-cohort study when
##the event is rare
ccpw <- function(alpha=0.95,hr,p1,N,pfail,q){
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
  pwr <- pnorm(
    Zalpha+sqrt(n)*theta*
      sqrt((p1*p2*pfail)/
          (q+(1-q)*pfail))
  )
  return(pwr)
}
