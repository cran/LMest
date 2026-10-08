forward_Multinom <- function(Pc,Piv,PI,n,k,Tv,fort=TRUE){

  TT = max(Tv)
  if(fort){
    out = .Fortran("forward_multinom",as.integer(sum(Tv)),as.integer(k),as.integer(n),Pc,Piv,
                   as.integer(Tv),as.integer(max(Tv)),PI,L = array(0,c(n,k,TT)),lkv = rep(0,n))
    L = out$L; lkv = out$lkv
  }else{
    TT = max(Tv)
    L = array(0,c(n,k,TT)); lkv = rep(0,n)
    j = 0
    for(i in 1:n){
      j = j+1
      L[i,,1] = Pc[j,]*Piv[i,]
      tmp = sum(L[i,,1]); lkv[i] = log(tmp)
      L[i,,1] = L[i,,1]/tmp
      for(t in 2:Tv[i]){
        j = j+1
        L[i,,t] = Pc[j,]*(L[i,,t-1]%*%PI[,,i,t])
        tmp = sum(L[i,,t]); lkv[i] = lkv[i] + log(tmp)
        L[i,,t] = L[i,,t]/tmp
      }
    }
  }
  out = list(L=L,lkv=lkv)
  return(out)

}
