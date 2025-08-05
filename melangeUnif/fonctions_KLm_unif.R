# divergence = "KLm"

# Mesure empirique du critère assosicé à KLm
Pn.m = function(X,pi,theta) mean(sapply(sapply(X,pi=pi,theta=theta,FUN=g),FUN=log))

# Première dérivée partielle p\r à pi de m
partial_pi.m = function(x,pi,theta) (df2(x,theta)-df1(x))/g(x,pi,theta)

bn = function(X,t1,t2){
  hatpi.1 = hat.pi(X,t1)
  hatpi.2 = hat.pi(X,t2)
  num = mean(sapply(X, pi = hatpi.1, theta = t1, FUN = partial_pi.m)*sapply(X, pi = hatpi.2, theta = t2, FUN = partial_pi.m))
  denom.1 = mean((sapply(X, pi = hatpi.1, theta = t1, FUN = partial_pi.m))**2)
  denom.2 = mean((sapply(X, pi = hatpi.2, theta = t2, FUN = partial_pi.m))**2)
  return(num/(denom.1*denom.2))
}

# 
temp.fun = function(x,pi,theta) (df2(x,theta) - df1(x))**2/(g(x,pi,theta)**2)
# 
stat.T = function(X,theta){
  
  hat.pi = function(theta) 
  
  hat.pi = optimize(f = Pn.m, X = X, theta = theta, interval = c(0,1), maximum = TRUE)$maximum
  an = 1/(mean(sapply(X, pi = hat.pi, theta = theta, FUN = temp.fun)))
  return(c(hat.pi, an, sqrt(n/an)*hat.pi))
}

# hat.pi
hatpi.function = function(matrixSample,eta){
  sample = c()
  for(i in c(1:nrow(matrixSample))){
    ligne = matrixSample[i,]
    n.moins = sum(ligne < eta)
    n.plus = sum(ligne > eta)
    hatpi = (n.moins - eta/(1-eta)*n.plus)/n
    sample = c(sample, hatpi)
  }
  # h.breaks = (1+eta/(1-eta))/n
  # assign("h.breaks", h.breaks, envir = .GlobalEnv)
  return(sample)
}

an.function = function(X, theta){
  n.moins = sum(ligne < eta)
  n.plus = sum(ligne > eta)
  hatpi = (n.moins - eta/(1-eta)*n.plus)/n
  sampleHATPI = c(sampleHATPI, hatpi)
  h.breaks = (1+eta/(1-eta))/n
  assign("h.breaks", h.breaks, envir = .GlobalEnv)
  # AN
  denom.1 = n.moins / (eta/(1-eta)+hatpi)**2
  denom.2 = n.plus / (1-hatpi)**2
  an = n / (denom.1 + denom.2)
}


# an
an = function(X,theta) 1 / (mean(sapply(X, pi = hat.pi(X,theta), theta = theta, FUN = temp.fun)))

rho = function(X,t1,t2) bn(X,t1,t2)/sqrt(an(X,t1)*an(X,t2))