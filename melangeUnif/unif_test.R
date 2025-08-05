rm(list = ls())
library(MASS)
library(pracma)
# setwd(dir = '/Users/cyril/cthommeret-phd/simulations/Rscripts/script_01')

source('specification.R')

# ETAPE 1: Spécification du mélange par la fonction spec();
# ETAPE 2: Définir D.Theta_2 # ETAPE 3: Définir Dtilde.Theta_2;

spec(distrib1.name = "unif", distrib2.name = "unif", distrib1.parameter = c(0,1), distrib2.parameter = c(0,"max"),indice.pos = 2)



# ETAPE 4: Faire tourner l'algorithme par la fonction algorithm_1.

stat.T = function(echantillon, eta){
  n = length(echantillon)
  n.moins = sum(echantillon < eta)
  n.plus = n - n.moins
  
  hatpi = (1+eta/(1-eta))*n.moins/n - eta/(1-eta)
  
  denom.1 = n.moins/(eta/(1-eta) + hatpi)**2
  denom.2 = n.plus/(1-hatpi)**2
  
  an = n/(denom.1+denom.2)
  
  t = sqrt(n/an)*hatpi
  
  return(t)
}


algorithm = function(divergence, nexp, n, n.tilde, I){
  
  assign("g", function(x,pi,theta) (1-pi)*df1(x) + pi*df2(x,theta), envir = .GlobalEnv)
  assign("r", rep(0,I), envir = .GlobalEnv)
  assign("p", runif(I), envir = .GlobalEnv)
  assign("n", n, envir = .GlobalEnv)
  assign("nexp", nexp, envir = .GlobalEnv)
  assign("n.tilde", n.tilde, envir = .GlobalEnv)
  assign("I", I, envir = .GlobalEnv)
  assign("nb2", length(Dtilde.Theta_2), envir = .GlobalEnv)
  
  if(divergence == "KLm"){
    #source('fonctions.KLm_unif.R')
    algorithm.KLm()
  }
  if(divergence == "chi2"){
    source('fonctions.chi2_unif.R')
    algorithm.chi2()
  }

}
