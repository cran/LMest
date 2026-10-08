comp_PI <- function(k,n,Tv,XX,GG,eta,fort=FALSE){

  nc = ncol(XX)-1
  TT = max(Tv)
  if(fort){
    out = .Fortran("comp_PI",as.integer(k),as.integer(n),as.integer(Tv),as.integer(TT),
                   as.integer(nc),nT=as.integer(sum(Tv)),XX=XX,GG,eta,PI=array(0,c(k,k,n,TT)))
    PI = out$PI
  }else{
    knc = (k-1)*(1+nc)
    PI = array(0,c(k,k,n,TT))
    j = 0
    for(i in 1:n){
      j = j+1
      for(t in 2:Tv[i]){
        j = j+1
        Xi = diag(k-1)%x%t(XX[j,])
        for(u in 1:k){
          ind = (u-1)*knc+(1:knc)
          tmp = c(GG[,,u]%*%Xi%*%eta[ind])
          tmp = exp(tmp-max(tmp))
          PI[u,,i,t] = tmp/sum(tmp)
        }
      }
    }
  }
  return(PI)

}
