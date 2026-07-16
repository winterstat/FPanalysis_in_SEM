converter <- function(per){
  
  res = per * 1*10^6
  
  return(res/100)
  
}

converter(per = .455)

# percentage increase
per_in <- function(N_v, O_v) {
  
  res <- (N_v - O_v)/O_v
  
  res * 100
  
}

per_in(N_v = 2.967, O_v = 1.115)

c <- 76187445