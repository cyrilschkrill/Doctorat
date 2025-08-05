# Melange Uniforme && KLm: Fonctions #
specificationDuMelange = function(nomDistribution.1, 
                                  nomDistribution.2, 
                                  parametrage.1, 
                                  nombreParametres.distribution_2 = 1,
                                  positionParam.2_inconnu = 1,
                                  valeurParam.2_connue = NaN){
  distributionList = c("exp", "norm", "lnorm", "weibull","unif")
  test = (nomDistribution.1 %in% distributionList) & (nomDistribution.1 %in% distributionList)
  if( test == FALSE ){ print("Erreur de saisie dans le nom d'une des composantes")}
  
  # PREMIERE COMPOSANTE # 
  if(length(parametrage.1) == 1){
    rf1 = function(n) get(paste0("r", nomDistribution.1))(n, parametrage.1) 
    df1 = function(x) get(paste0("d", nomDistribution.1))(x, parametrage.1)
  }
  if(length(parametrage.1) == 2){
    rf1 = function(n) get(paste0("r", nomDistribution.1))(n, parametrage.1[1], parametrage.1[2])
    df1 = function(x) get(paste0("d", nomDistribution.1))(x, parametrage.1[1], parametrage.1[2])
  }
  assign("rf1", rf1, envir = .GlobalEnv)
  assign("df1", df1, envir = .GlobalEnv)
  
  # DEUXIEME COMPOSANTE #
  if(nombreParametres.distribution_2 == 1){
    rf2 = function(n,eta) get(paste0("r", nomDistribution.2))(n, eta)
    df2 = function(x,eta) get(paste0("d", nomDistribution.2))(x, eta)
  }
  else{ # nombreParametres.distribution_2 == 2
    if(nombreParametres.distribution_2 != 2) return("Erreur: Le nombre de paramètres de la composante ne peut que prendre les valeurs 1 ou 2.")
    if(positionParam.2_inconnu == 1){
      rf2 = function(n,eta) get(paste0("r", nomDistribution.2))(n, eta, valeurParam.2_connue)
      df2 = function(x,eta) get(paste0("d", nomDistribution.2))(x, eta, valeurParam.2_connue)
    }
    else{ # positionParam.2_inconnu == 2
      if(positionParam.2_inconnu != 2) return("Erreur: La position du paramètre inconnu de la seconde composante.")
      rf2 = function(n,eta) get(paste0("r", nomDistribution.2))(n, valeurParam.2_connue, eta)
      df2 = function(x,eta) get(paste0("d", nomDistribution.2))(x, valeurParam.2_connue, eta)
    }
  }
  assign("rf2", rf2, envir = .GlobalEnv)
  assign("df2", df2, envir = .GlobalEnv)
  
  # FONCTIONS DU MELANGE #
  rmix = function(n,pi,eta){
    output = sample(x = c(0,1), size = n, replace = TRUE, prob = c(1-pi,pi))
    nb1 = sum(output)
    output[output == 1] = rf2(nb1,eta)
    output[output == 0] = rf1(n-nb1)
    return(output)
  }
  dmix = function(x,pi,eta) (1-pi)*df1(x) + pi*df2(x,eta)
  assign("rmix", rmix, envir = .GlobalEnv)
  assign("dmix", dmix, envir = .GlobalEnv)
}

MATRICE_ECHANTILLONNAGE = function(nexp, n, pi.star, eta.star){
  assign('nexp',nexp,envir = .GlobalEnv)
  if(!exists("rmix")) "La spécification du mélange n'a pas opérée."
  
  matrice.dechantillonnage = matrix(data = rmix(n = n*nexp, pi = pi.star, eta = eta.star), nrow = nexp, byrow = TRUE)
  assign("matrice.dechantillonnage", matrice.dechantillonnage, envir = .GlobalEnv)
  
  # return(matrice.dechantillonnage)
}

CALCUL_STATISTIQUE = function(echantillon, eta){
  n = length(echantillon)
  n.moins = sum(echantillon < eta)
  hatpi = (1+eta/(1-eta))*n.moins/n - eta/(1-eta)
  
  n.plus = n - n.moins
  
  if(n.plus == 0 | n.moins == 0){
    an = 1/(1-eta)**2
  }
  else{
    denom.1 = n.moins/(eta/(1-eta) + hatpi)**2
    denom.2 = n.plus/(1-hatpi)**2
    an = n/(denom.1 + denom.2)
  }
  
  s = sqrt(n/an)*hatpi
  
#  return(data.frame(hatpi,an,s))
  return(c(hatpi,an,s))
}

MATRICE_STATISTIQUE = function(matrice.dechantillonnage, eta){
  matrice.statistique = t(apply(X = matrice.dechantillonnage, eta = eta, FUN = CALCUL_STATISTIQUE, MARGIN = 1))
  
  matrice.statistique = data.frame(matrice.statistique)
  names(matrice.statistique) = c("hatpi","an","s")
  
  assign("matrice.statistique", matrice.statistique, envir = .GlobalEnv)
  #return(matrice.statistique)
}

VECTEUR_T = function(matrice.dechantillonnage, discretisation.Theta_2){
  vecteur.s = sapply(X = discretisation.Theta_2, FUN = (function(eta) MATRICE_STATISTIQUE(matrice.dechantillonnage, eta)$s) )
  assign("vecteur.s", vecteur.s, envir = .GlobalEnv)
  vecteur.t = apply(X = vecteur.s, FUN = max, MARGIN = 1)
  assign("vecteur.t", vecteur.t, envir = .GlobalEnv)
}

