###########
# ETAPE 0 #
###########

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

###########
# ETAPE 1 #
###########

MATRICE_ECHANTILLONNAGE = function(nexp, n, pi.star, eta.star){
  assign('nexp',nexp,envir = .GlobalEnv)
  if(!exists("rmix")) "La spécification du mélange n'a pas opérée."
  
  matrice.dechantillonnage = matrix(data = rmix(n = n*nexp, pi = pi.star, eta = eta.star), nrow = nexp, byrow = TRUE)
  assign("matrice.dechantillonnage", matrice.dechantillonnage, envir = .GlobalEnv)
  
  # return(matrice.dechantillonnage)
}

###########
# ETAPE 2 #
###########

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

###########
# ETAPE 3 #
###########

MATRICE.TKK.TILDE = function(n.tilde, matrice.dechantillonnage, discretisation.Theta_2.tilde){
  
  nb2 = length(discretisation.Theta_2.tilde)
  vecteurDeMatrice = apply(X = matrice.dechantillonnage, discretisation.Theta_2.tilde = discretisation.Theta_2.tilde, FUN = matriceDeVarianceCovariance,  MARGIN = 1)
  
  output = c()
  for(matriceDeCovariance in vecteurDeMatrice){
    echantillon.ttilde = (function(SIGMA) replicate(n.tilde, max.Y(SIGMA,nb2)))(matriceDeCovariance)
    output = rbind(output, echantillon.ttilde)
  }
  assign("matrice.tkk_tilde",output, envir = .GlobalEnv)
}

max.Y = function(matriceCovariance,nb2){
  Y = mvrnorm(n = 1, mu = rep(x = 0, nb2), Sigma = matriceCovariance)
  t = max(Y)
  return(t)
}

matriceDeVarianceCovariance = function(echantillon,discretisation.Theta_2.tilde){
  nb2 = length(discretisation.Theta_2.tilde)
  output = c()
  for(i in c(1:nb2)){
    ligne = c()
    for(j in c(1:nb2)){
      if( j < i){
        ligne = c(ligne,0)
      }
      else{
        t1 = discretisation.Theta_2.tilde[i]
        t2 = discretisation.Theta_2.tilde[j]
        iter = rho(echantillon,t1,t2)
        ligne = c(ligne,iter)
      }
    }
    output = rbind(output,ligne)
  }
  temp.mat = t(output)
  diag(temp.mat) = rep(0,nb2)
  output = output + temp.mat
  print(paste0(compteur," /",nexp))
  assign("compteur",compteur+1,envir = .GlobalEnv)
  return(as.data.frame(output))
}


# Première dérivée partielle p\r à pi de m
partial_pi.m = function(x,pi,eta){
  if(x < eta | x == eta){
    output = (1 - eta)/(eta + pi*(1 - eta))
  }
  else{
    output = 1/(pi - 1)
  }
  return(output)
}

HATPI = function(echantillon,eta){
  n.moins = sum(echantillon < eta)
  hatpi = (1+eta/(1-eta))*n.moins/length(echantillon) - eta/(1-eta)
  return(hatpi)
}

rho = function(echantillon,eta.1,eta.2){
  n = length(echantillon)
  if(eta.1 == eta.2){
    hatpi = HATPI(echantillon,eta.1)
    bn = 1/mean((sapply(echantillon, pi = hatpi, eta = eta.1, FUN = partial_pi.m))**2)
    
    n.moins = sum(echantillon < eta.1)
    n.plus = n - n.moins
    if(n.plus == 0 | n.moins == 0){
      an = 1/(1-eta.1)**2
    }
    else{
      denom.1 = n.moins/(eta.1/(1-eta.1) + hatpi)**2
      denom.2 = n.plus/(1-hatpi)**2
      an = n/(denom.1 + denom.2)
    }
    return(bn/an) 
  }
  else{
    hatpi.1 = HATPI(echantillon,eta.1)
    hatpi.2 = HATPI(echantillon,eta.2)
    
    bn = mean(sapply(echantillon, pi = hatpi.1, eta = eta.1, FUN = partial_pi.m)*sapply(echantillon, pi = hatpi.2, eta = eta.2, FUN = partial_pi.m)) / ( mean((sapply(echantillon, pi = hatpi.1, eta = eta.1, FUN = partial_pi.m))**2)*mean((sapply(echantillon, pi = hatpi.2, eta = eta.2, FUN = partial_pi.m))**2))
    
    # AN _1
    n.moins = sum(echantillon < eta.1)
    n.plus = n - n.moins
    if(n.plus == 0 | n.moins == 0){
      an.1 = 1/(1-eta.1)**2
    }
    else{
      denom.1 = n.moins/(eta.1/(1-eta.1) + hatpi.1)**2
      denom.2 = n.plus/(1-hatpi.1)**2
      an.1 = n/(denom.1 + denom.2)
    }
    
    # AN _2
    n.moins = sum(echantillon < eta.2)
    n.plus = n - n.moins
    if(n.plus == 0 | n.moins == 0){
      an.2 = 1/(1-eta.2)**2
    }
    else{
      denom.1 = n.moins/(eta.2/(1-eta.2) + hatpi.2)**2
      denom.2 = n.plus/(1-hatpi.2)**2
      an.2 = n/(denom.1 + denom.2)
    }
    return(bn/sqrt(an.1*an.2))
  }
  
}

###########
# ETAPE 4 #
###########

VECTEUR_REJET = function(vecteur.t, matrice.tkk_tilde,p){
  q=1-p
  I = length(p)
  nb.rejet = rep(0,I)
  n.tilde = ncol(matrice.tkk_tilde)
  nexp = length(vecteur.t)
  
  for(i in c(1:I)){
    for(k in c(1:nexp)){
      ordered.tkk_tilde = matrice.tkk_tilde[k,][order(matrice.tkk_tilde[k,])]
      nb.rejet[i] = nb.rejet[i] + 1*(vecteur.t[k] > ordered.tkk_tilde[as.integer(n.tilde*q[i])+1] | vecteur.t[k] == ordered.tkk_tilde[as.integer(n.tilde*q[i])+1])
    }
  }
  return(nb.rejet/nexp)
}


