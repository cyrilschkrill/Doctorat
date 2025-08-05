stat.t = function(echantillon, eta){
  n = length(echantillon)
  n.moins = sum(echantillon < eta)
  n.plus = n - n.moins
  # print("n.moins")
  # print(n.moins)
  hatpi = (1+eta/(1-eta))*n.moins/n - eta/(1-eta)
  # print("hatpi")
  # print(hatpi)
  denom.1 = n.moins/(eta/(1-eta) + hatpi)**2
  denom.2 = n.plus/(1-hatpi)**2
  
  an = n/(denom.1+denom.2)
  # print("an")
  # print(an)
  
  t = sqrt(n/an)*hatpi
  
  return(t)
}
t = function(echantillon){
  l = sapply(X = discretisationTheta, FUN = stat.t, echantillon = echantillon)
  l = l[!is.na(l)] #supprime les incohérences aux bords (n.moins = 0 et n.plus = 1)
  return(max(l))
}

algorithm.KLm = function(n,nexp){
  matrice.dechantillonnage = matrix(data = rf1(n*nexp), byrow = TRUE, nrow = nexp)
  vecteur.tk = apply(X = matrice.dechantillonnage, FUN = t, MARGIN = 1)
  
  return(list("vecteur.tk" = vecteur.tk, "MatCov" = MatCov))
}

spec(distrib1.name = "unif", distrib2.name = "unif", distrib1.parameter = c(0,1), distrib2.parameter = list(known.parameter=0,unknown.parameter = "max"), indice.pos = 2)



discretisationTheta = runif(100)
algorithm.KLm(2,5)




algorithm.KLm = function(){
  vecteur.tk = c()
  for(k in c(1:nexp)){
    tk = max(sapply(X = discretisationTheta,
                    FUN = stat.t,
                    echantillon = rmix(n)))
    vecteur.tk = c(vecteur.tk,tk)
    
  }
}


algorithm_1.KLm = function(){
  
  vect.tk = c()
  mat.T = c()
  for(k in c(1:nexp)){
    # Etape 1
    X = rf1(n)
    Tk = sapply(D.Theta_2, X = X, FUN = stat.T)#################################
    # Etape 2
    tk = max(Tk[3,])
    vect.tk = c(vect.tk,tk)
    
    # Matrice de covariance
    SIG = c()
    for(i in c(1:nb2)){
      ligne = c()
      for(j in c(1:nb2)){
        if( j < i){
          ligne = c(ligne,0)
        }
        else{
          t1 = Dtilde.Theta_2[i]
          t2 = Dtilde.Theta_2[j]
          iter = rho(X,t1,t2) ##################################################
          ligne = c(ligne,iter)
        }
      }
      SIG = rbind(SIG,ligne)
    }
    temp.mat = t(SIG)
    diag(temp.mat) = rep(0,nb2)
    SIG = SIG + temp.mat
    assign("matCov", SIG, envir = .GlobalEnv)
    
    # Etape 3
    ligne = c()
    for(ktilde in c(1:n.tilde)){
      # (a)
      Y = mvrnorm(n = nb2, mu = rep(0,nb2), Sigma = SIG)
      # (b)
      tktilde = max(Y)
      ligne = c(ligne,tktilde)
    }
    mat.T = rbind(mat.T,ligne)
    
    # Etape 4
    for(i in c(1:I)){
      if(tk < ligne[as.integer(n.tilde*p[i])+1]){
        r[i] = r[i]
      }
      else{ 
        r[i] = r[i] + 1/nexp
      }
    }
  }
  assign("vect.tk", vect.tk, envir = .GlobalEnv)
  assign("mat.T", mat.T, envir = .GlobalEnv)
  return(c(vect.tk,mat.T,r))
}
