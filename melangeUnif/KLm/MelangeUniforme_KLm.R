library(pracma)
library(MASS)
# Melange Uniformes dans le cadre KLm
temp.STATISTIQUE = function(echantillon, eta){
  n.moins = sum(echantillon < eta)
  hatpi = (1+eta/(1-eta))*n.moins/length(echantillon) - eta/(1-eta)
  
  n.plus = n - n.moins

  if(n.plus == 0 | n.moins == 0){
    an = 1/(1-eta)**2
  }
  else{
    denom.1 = n.moins/(eta/(1-eta) + hatpi)**2
    denom.2 = n.plus/(1-hatpi)**2
    an = n/(denom.1 + denom.2)
  }
  t = sqrt(n/an)*hatpi
  
  return(c(hatpi,an,t))
}

STATISTIQUE = function(matrice.dechantillonnage, eta){
  output = t(apply(X = matrice.dechantillonnage, MARGIN = 1, FUN = temp.STATISTIQUE, eta = eta))
  output = data.frame(output)
  names(output) = c("hatpi","an","t")
  return(output)
}

VECTEUR.TK = function(matrice.dechantillonnage){
  temp = function(eta) STATISTIQUE(matrice.dechantillonnage, eta)$t
  t = sapply(X = discretisation.Theta_2, FUN = temp)
  output = apply(X = t, MARGIN = 1,FUN = max)
  return(t(output))
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

BN = function(X,t1,t2){
  hatpi.1 = HATPI(X,t1)
  hatpi.2 = HATPI(X,t2)
  
  num = mean(sapply(X, pi = hatpi.1, eta = t1, FUN = partial_pi.m)*sapply(X, pi = hatpi.2, eta = t2, FUN = partial_pi.m))
  
  denom.1 = mean((sapply(X, pi = hatpi.1, eta = t1, FUN = partial_pi.m))**2)
  denom.2 = mean((sapply(X, pi = hatpi.2, eta = t2, FUN = partial_pi.m))**2)
  
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
  return(as.data.frame(output))
}

vecteurLigne = function(matriceCovariance){
  ligne = replicate(n = n.tilde, max.Y(matriceCovariance))
  return(ligne)
}

max.Y = function(matriceCovariance){
  Y = mvrnorm(n = nb2, mu = rep(x = 0, nb2), Sigma = matriceCovariance)
  t = max(Y)
  return(t)
}

MATRICE.TKK.TILDE = function(matrice.dechantillonnage){
  vecteurDeMatrice = apply(X = matrice.dechantillonnage, MARGIN = 1, FUN = matriceDeVarianceCovariance)
  output = c()
  for(matrice in vecteurDeMatrice){
    ligne = vecteurLigne(matrice)
    output = rbind(output, ligne)
  }
  return(output)
}

VECTEUR.REJET = function(vecteurTk, matriceTkk.tilde){
  r = rep(0,I)
  for(k in c(1:length(vecteurTk))){
    for(i in c(1:I)){
      ligne = matriceTkk.tilde[k,]
      if(vecteurTk[k] < ligne[as.integer(n.tilde*p[i])+1]){
        r[i] = r[i]
      }
      else{ 
        r[i] = r[i] + 1/nexp
      }
    }
  }
  return(r)
}

# EXECUTE PROCESS
n = 30
nexp = 100
matrice.dechantillonnage = matrix(data = runif(n = n*nexp), nrow = nexp, byrow = TRUE)

n.tilde = 100

discretisation.Theta_2 = linspace(0.01,0.999,30)

discretisation.Theta_2.tilde = discretisation.Theta_2
nb2 = length(discretisation.Theta_2.tilde)

I=20
p = runif(I)

#
tk = VECTEUR.TK(matrice.dechantillonnage)
matrice.tkk.tilde = MATRICE.TKK.TILDE(matrice.dechantillonnage)
r = VECTEUR.REJET(tk, matrice.tkk.tilde)

# print(tk)
# print(matrice.tkk.tilde)
print(r)



























    