# setwd(dir = '~/Desktop/phd/melangeUnif/chi2/scripts')
rm(list = ls())
source('functionScript_chi2.R')
assign("compteur",1,envir = .GlobalEnv)

specificationDuMelange(nomDistribution.1 = "unif", nomDistribution.2 = "unif", parametrage.1 = c(0,1), nombreParametres.distribution_2 = 2, positionParam.2_inconnu = 2, valeurParam.2_connue = 0)

eta.star = 0.33
pi.star = 0.5

# 1
MATRICE_ECHANTILLONNAGE(nexp = 10, n = 15, pi.star = 0.5, eta.star = 0.3)

# 2
library(pracma)
discretisation.Theta_2 = linspace(0.1,0.9,20)
VECTEUR_T(matrice.dechantillonnage = matrice.dechantillonnage, discretisation.Theta_2 = discretisation.Theta_2)

# 3
library(MASS)
discretisation.Theta_2.tilde = discretisation.Theta_2
MATRICE.TKK.TILDE(n.tilde = 7, matrice.dechantillonnage = matrice.dechantillonnage, discretisation.Theta_2.tilde = discretisation.Theta_2.tilde)

# 4
p = c(0.01,0.05,0.1,0.2,0.5)
VECTEUR_REJET(vecteur.t = vecteur.t, matrice.tkk_tilde = matrice.tkk_tilde, p = p)