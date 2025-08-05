# Melange Uniformes dans le cadre KLm
n = 30
nexp = 5
matrice.dechantillons = matrix(data = runif(n = n*nexp), nrow = nexp, byrow = TRUE)

STATISTIQUE = function(echantillon, eta){
  n.moins = sum(echantillon < eta)
  hatpi = (1+eta/(1-eta))*n.moins/length(echantillon) - eta/(1-eta)
  
  n.plus = n - n.moins
  denom.1 = n.moins/(eta/(1-eta) + hatpi)**2
  denom.2 = n.plus/(1-hatpi)**2
  
  an = n/(denom.1 + denom.2)
  
  t = sqrt(n/an)*hatpi
  
  return(data.frame(hatpi,an,t))
}