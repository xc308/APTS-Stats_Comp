#======================#
# Numerical Integration
#======================#


# find ways to approximate an intergral of a function 
# with minimum number of function evaluations


#---------------------------------#
# Quadrature Rules: Mid-point rule
#---------------------------------#

# integration wrt 2 or up to 4 variables, either use
  # differential eq solo methods or 
  # quadrature rules
    # evaluate the function (integrand) at some set of points
    # in the integration domain & combine resulting values
    # to approximate the integral


# e.g int from 0 to 1 of exp(x)
I.true <- exp(1) - exp(0)

N <- 10
1:N
I.mp <- sum(exp(0 + (1:N - .5) / N)) / N

c(I.mp, I.true, (I.true - I.mp) / I.true) 
# [1] 1.7175660865 1.7182818285 0.0004165452

# .4% error not bad
# but want higher accuracy

N1 <- 100
N2 <- 1000

Error_band <- function(N) {
  I.mp <- sum(exp((1:N - .5) / N)) / N
  return((I.true - I.mp) / I.true)
}

Error_band(N1)  # [1] 4.166655e-06
Error_band(N2)  # [1] 4.166667e-08


## so each 10 fold increase in the number of evaluations 
# buys us a 100 fold reduction in the approximation error. 
# and this error is not a random error, but a real bias. 


#------------------------------------#
# Quadrature Rules:  Guass Quadrature
#------------------------------------#

# statmod package has gauss.quad function producing
# weights wi and nodes xi for a fixed weight function w(x)
# st. int from a to b of g(x)w(x) = sum from i = 1 to n wig(xi)

# it assumes a= -1, b = 1, but can re-scale to actual interval

install.packages("statmod")
library(statmod)

N <- 10
gq <- gauss.quad(N)

gq

sum(gq$weights) #2


I.gq <- sum(gq$weights * exp((gq$nodes + 1)/2))/2

(I.true - I.gq) / I.true  #[1] -1.292248e-15


Err_gq <- function(N) {
  
  gq <- gauss.quad(N)
  I.gq <- sum(gq$weights * exp((gq$nodes + 1) / 2)) / 2
  return((I.true - I.gq) / I.true)
 
}

N1 <- 10
N2 <- 5
N3 <- 3

Err_gq(N1)  # [1] -1.292248e-15
Err_gq(N2)  # [1] 3.80567e-13
Err_gq(N3)  # [1] 4.795992e-07

## even just 3 nodes is more accurate than mid-point rule


































































