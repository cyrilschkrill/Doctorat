Pm = function(echantillon, pi, eta){
  a = 1 + pi*(1/eta-1)
  b = 1 - pi
  
  n.moins = sum(echantillon<eta)
  n = length(echantillon)
  n.plus = n - n.moins
  
  A = (b*eta + a*(1-eta))/(a*b)
  B = 0.5*(1+(n.moins/a**2 + n.plus/b**2)/n)
    
  output = A - B
  return(output)
}
#################
# Test script 1 #
#################
#Pm(X = runif(100), pi = 0.5, eta = 0.3)
# curve(Pm(X = runif(1e3), pi = x, eta = 0.8),from = 0, to = 1)
# 
# specificationDuMelange(nomDistribution.1 = "unif", nomDistribution.2 = "unif", parametrage.1 = c(0,1), nombreParametres.distribution_2 = 2, positionParam.2_inconnu = 2, valeurParam.2_connue = 0)
# 
# X.test = rmix(1e3,0.5,0.5)
# hist(X.test, breaks = 25)
# curve(Pm(X = X.test, pi = x, eta = 0.9),from = 0, to = 0.8)

HATPI = function(echantillon, eta){
  temp.fun = function(pi) Pm(echantillon = echantillon, pi = pi, eta = eta)
  return(optimise(f = temp.fun, interval = c(0,1), maximum = TRUE)$maximum)
}  

library(pracma)
HATPI.2 = function(echantillon, eta){
  temp.fun = function(pi) Pm(echantillon = echantillon, pi = pi, eta = eta)
  newtonRaphson(fun = temp.fun, x0 = 0.5)
}  
