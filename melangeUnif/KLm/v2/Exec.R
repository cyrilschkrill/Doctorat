setwd(dir = '/Users/cyril/cthommeret-phd/simulations/melangeUnif/KLm/v2')
rm(list = ls())
source('Fonctions.R')
assign("compteur",1,envir = .GlobalEnv)

specificationDuMelange(nomDistribution.1 = "unif", nomDistribution.2 = "unif", parametrage.1 = c(0,1), nombreParametres.distribution_2 = 2, positionParam.2_inconnu = 2, valeurParam.2_connue = 0)

x = rmix(1e3,pi = 0.5, eta = 0.33)
hist(x,breaks = 100)

# CALCUL_STATISTIQUE(x, eta = 0.3)

MATRICE_ECHANTILLONNAGE(nexp = 10, n = 1e3, pi.star = 0, eta.star = 0.3)
MATRICE_STATISTIQUE(matrice.dechantillonnage, eta = 0.5)
hist(matrice.statistique$hatpi,breaks=5)

library(pracma)
discretisation.Theta_2 = linspace(0.001,0.9,50)
discretisation.Theta_2.tilde = discretisation.Theta_2

VECTEUR_T(matrice.dechantillonnage,discretisation.Theta_2)
hist(vecteur.t,breaks=100,prob=TRUE)
# OK tout baigne !

source(file = 'ttilde.R')
library(MASS)
MATRICE.TKK.TILDE(n.tilde = 50, matrice.dechantillonnage,discretisation.Theta_2.tilde)

source(file = 'rejet2.R')
VECTEUR_REJET(vecteur.t, matrice.tkk_tilde,p=c(0.05,0.10,0.20,0.30))

  

