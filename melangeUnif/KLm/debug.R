discretisation.Theta_2.tilde = discretisation.Theta_2
nb2 = length(discretisation.Theta_2.tilde)


partial_pi.m = function(x,pi,theta){
  if(x < eta | x == eta){
    output = (1 - eta)/(eta + pi*(1 - eta))
  }
  else{
    output = 1/(pi - 1)
  }
  return(output)
}

BN = function(X,t1,t2){
  hatpi.1 = HATPI(X,t1)
  hatpi.2 = HATPI(X,t2)
  num = mean(sapply(X, pi = hatpi.1, theta = t1, FUN = partial_pi.m)*sapply(X, pi = hatpi.2, theta = t2, FUN = partial_pi.m))
  denom.1 = mean((sapply(X, pi = hatpi.1, theta = t1, FUN = partial_pi.m))**2)
  denom.2 = mean((sapply(X, pi = hatpi.2, theta = t2, FUN = partial_pi.m))**2)
  return(num/(denom.1*denom.2))
}

HATPI = function(echantillon,eta){
  n.moins = sum(echantillon < eta)
  hatpi = (1+eta/(1-eta))*n.moins/length(echantillon) - eta/(1-eta)
  return(hatpi)
}

AN = function(echantillon,eta){
  hatpi = HATPI(echantillon,eta)
  
  n.moins = sum(echantillon < eta)
  n.plus = n - n.moins
  
  if(n.plus == 0 | n.moins == 0){
    an = 1/(1-eta)**2
  }
  else{
    denom.1 = n.moins/(eta/(1-eta) + hatpi)**2
    denom.2 = n.plus/(1-hatpi)**2
    an = n/(denom.1 + denom.2)
  }
  
  return(an)
}

rho = function(X,t1,t2) BN(X,t1,t2)/sqrt(AN(X,t1)*AN(X,t2))

matriceDeVarianceCovariance = function(echantillon){
  
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
  return(output)
}


# X = matrice.dechantillonnage[1,]
# 
# eta = 0.02
# # hatpi = -0.02040816
# 
# 
# HATPI = function(echantillon,eta){
#   n.moins = sum(echantillon < eta)
#   hatpi = (1+eta/(1-eta))*n.moins/length(echantillon) - eta/(1-eta)
#   return(hatpi)
# }
# 
# 
# 
# AN = function(echantillon,eta){
#   hatpi = HATPI(echantillon,eta)
#   
#   n.moins = sum(echantillon < eta)
#   n.plus = n - n.moins
#   denom.1 = n.moins/(eta/(1-eta) + hatpi)**2
#   denom.2 = n.plus/(1-hatpi)**2
#   
#   an = n/(denom.1 + denom.2)
#   
#   return(data.frame(hatpi,n.plus,n.moins,denom.1,denom.2))
# }
# AN(X,0.02)

