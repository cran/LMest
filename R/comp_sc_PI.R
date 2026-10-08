comp_sc_PI <- function(k,n,Tv,XX,GG,W2,eta,fort=FALSE){

  nc = ncol(XX)-1
  knc = (k-1)*(1+nc)
  if(fort){
    nc = ncol(XX)-1
    out = .Fortran("comp_sc_PI",as.integer(k),as.integer(n),as.integer(Tv),as.integer(max(Tv)),
                   as.integer(nc),nT=as.integer(sum(Tv)),XX=XX,GG,W2,eta,sv2 = rep(0,k*knc),
                   J2 = matrix(0,k*knc,k*knc))
    sv2 = out$sv2; J2 = out$J2
  }else{
    sv2 = rep(0,k*knc)
    J2 = matrix(0,k*knc,k*knc)
    j = 0
    for(i in 1:n){
      j = j+1
      for(t in 2:Tv[i]){
        j = j+1
        Xit = diag(k-1)%x%t(XX[j,])
        for(u in 1:k){
          ind = (u-1)*knc+(1:knc)
          GXit = GG[,,u]%*%Xit
          tmp = c(GXit%*%eta[ind])
          tmp = exp(tmp-max(tmp))
          tmp = tmp/sum(tmp)
          tot = sum(W2[u,,i,t])
          sv2[ind] = sv2[ind]+t(GXit)%*%(W2[u,,i,t]-tot*tmp)
          Omi = diag(tmp)-tmp%o%tmp
          J2[ind,ind] = J2[ind,ind]+tot*t(GXit)%*%Omi%*%GXit
        }
      }
    }
  }
  out = list(sv2=sv2,J2=J2)
  return(out)

}