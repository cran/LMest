backward_Multinom <- function(L,Pc,PI,n,k,Tv,fort=FALSE){

  TT = max(Tv)
  if(fort){
    out = .Fortran("backward_multinom",as.integer(n),as.integer(k),as.integer(Tv),as.integer(max(Tv)),
                   as.integer(sum(Tv)),L,Pc,PI,W=array(0,c(n,k,TT)),W2=array(0,c(k,k,n,TT)))
    W = out$W; W2 = out$W2
  }else{
    TT = max(Tv)
    M = array(0,c(n,k,TT))
    for(i in 1:n) M[,,Tv[i]] = 1
    W2 = array(0,c(k,k,n,TT))
    j = sum(Tv)
    for(i in n:1){
      for(t in (Tv[i]-1):1){
        j = j-1
        M[i,,t] = PI[,,i,t+1]%*%(Pc[j+1,]*M[i,,t+1])
        M[i,,t] = M[i,,t]/sum(M[i,,t])
        Tmp = (L[i,,t]%o%(Pc[j+1,]*M[i,,t+1]))*PI[,,i,t+1]
        W2[,,i,t+1] = Tmp/sum(Tmp)
      }
      j = j-1
    }
    W = L*M
    for(t in 1:TT){
      ind = which(Tv>=t)
      W[ind,,t] = (1/rowSums(W[ind,,t]))*W[ind,,t]
    }
  }
  out = list(W=W,W2=W2)
  return(out)
  
}