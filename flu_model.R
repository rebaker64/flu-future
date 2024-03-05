sirs_ode <- function(time, state ,theta) {
  #browser()
  ## Parameters:
  qList <- theta[["qList"]]
  D <- theta[["D"]]
  L <- theta[["L"]]
  R0max <- theta[["R0max"]]
  R0min <- theta[["R0min"]]
  ause <- theta[["ause"]]

  #derived
  quse <- qList[time]
  
  ## States:
  S <- state["S"]
  I <- state["I"]
  R <- state["R"]
  N <- S + I + R
  
  ## ODEs:
  R0 = exp(ause*quse + log(R0max - R0min)) + R0min 
  beta = R0/D
  dS <- (R/L) -beta * S * I/N 
  dI <- beta * S * I/N - (I/D) 
  dR <- (I/D) - (R/L) 
  
  return(list(c(dS, dI, dR)))
}

##### example parameters
#paras = list(D = 4/7, L = 40, R0min = 1.2, R0max = 3, qList = qList)


