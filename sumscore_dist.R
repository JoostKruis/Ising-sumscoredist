## Copyright 2019, Joost Kruis, All rights reserved.


# description -------------------------------------------------------------

Calculate the analytical sumscore distribution for a 0/1 Ising model. 

## INPUT

# For a system with n nodes:

# omega : the connectivity matrix [n x n]
# mu : external magnetic field [n] 
# beta : the inverse temperature [1]

# function ----------------------------------------------------------------

ising.sum.dist <- function(omega, mu, beta = 1)
{
  
  n = length(mu) # number of magnets
  
  x.all <- as.matrix(expand.grid(rep(list(as.integer(c(0,1))),n))) # all configurations
  
  boltzmann <- function(x,b=1){exp(-b*x)} # boltzmann function
  hamiltonian <- function(x,mu,omega){(1/2)*(-t(x)%*%omega%*%x) - t(x)%*%mu} # hamiltonian function
  
  E = apply(x.all,1,hamiltonian,mu=mu,omega=omega) # calculate energy for each configuration
  B = boltzmann(E,b = beta) # apply boltzmann function to energies
  Z <- sum(B)
  
  Y = cbind(B/Z,rowSums(x.all)) # combine probabilities with rpwsums

  colnames(Y)= c("P","X")
  
  Y %<>% data.frame %>% group_by(X) %>% summarize(P.sum = sum(P)) %>% data.frame
  
  print(round(Y,3))
  
  return(Y)
}


# example ----------------------------------------------------------------

# number of nodes
n = 5 

# connectivity matrix
omega <- matrix(1,n,n) 
diag(omega) = 0

# external magnetic field
mu <- rep(0,n)

# inverse temperature
beta = 1

# sum score distribution
res <- ising.sum.dist(omega, mu, beta)

# example output (print)

# X = sumscore
# P.sum = probability sumscore

##   X P.sum
## 1 0 0.000
## 2 1 0.000
## 3 2 0.001
## 4 3 0.008
## 5 4 0.083
## 6 5 0.907
