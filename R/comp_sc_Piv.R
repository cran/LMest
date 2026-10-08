comp_sc_Piv <- function(n,k,XX1,G,W,la,fort){
  
  nc1 = ncol(XX1)-1
  knc1 = (k-1)*(1+nc1)
  if(fort){
    out = .Fortran("comp_sc_Piv",as.integer(n),as.integer(k),as.integer(nc1),XX1,G,W[,,1],la,
                   sv1 = rep(0,knc1),J1 = matrix(0,knc1,knc1))
    sv1 = out$sv1; J1 = out$J1
  }else{
    sv1 = rep(0,knc1)
    J1 = matrix(0,knc1,knc1)
    for(i in 1:n){
      Xi = diag(k-1)%x%t(XX1[i,])
      GXi = G%*%Xi
      tmp = c(GXi%*%la)
      tmp = exp(tmp-max(tmp))
      tmp = tmp/sum(tmp)
      sv1 = sv1+t(GXi)%*%(W[i,,1]-tmp)
      Omi = diag(tmp)-tmp%o%tmp
      J1 = J1+t(GXi)%*%Omi%*%GXi
    }
  }
  out = list(sv1=sv1,J1=J1)
  return(out)

}
