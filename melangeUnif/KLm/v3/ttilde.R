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









rho = function(echantillon,eta.1,eta.2){
  n = length(echantillon)
  if(eta.1 == eta.2){
    hatpi = HATPI(echantillon,eta.1)
    an = AN(echantillon = echantillon, hatpi = hatpi, eta = eta)
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





