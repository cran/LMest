comp_Piv <- function(n,k,XX1,G,la,fort=FALSE){

  if(fort){
    nc1 = ncol(XX1)-1
    out = .Fortran("comp_Piv",as.integer(n),as.integer(k),as.integer(nc1),XX1,G,la,
                   Piv = matrix(0,n,k))
    Piv = out$Piv
  }else{
    Piv = matrix(0,n,k)
    for(i in 1:n){
      Xi = diag(k-1)%x%t(XX1[i,])
      tmp = c(G%*%Xi%*%la)
      tmp = exp(tmp-max(tmp))
      Piv[i,] = tmp/sum(tmp)
    }
  }
  # Piv = pmin(Piv,1-10^-100)
  # Piv = pmax(Piv,10^-100)
  return(Piv)

}