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