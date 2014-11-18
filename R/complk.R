complk <- function(S,yv,piv,Pi,Psi,k){

# Preliminaries
  	sS = dim(S)
  	ns = sS[1]
  	TT = sS[2]
  	if(length(sS)==2) r = 1 else r = sS[3]
# Compute log-likelihood
	Phi = array(1,c(ns,k,TT)); L = array(0,c(ns,k,TT))
	if(r==1) Phi[,,1] = Phi[,,1]*Psi[S[,1]+1,,1] else for(j in 1:r) Phi[,,1] = Phi[,,1]*Psi[S[,1,j]+1,,j]
  	L[,,1] = Phi[,,1]%*%diag(piv)
  	for(t in 2:TT){
  		if(r==1) Phi[,,t] = Phi[,,t]*Psi[S[,t]+1,,1] else for(j in 1:r) Phi[,,t] = Phi[,,t]*Psi[S[,t,j]+1,,j]
   		L[,,t] = Phi[,,t]*(L[,,t-1]%*%Pi[,,t])
  	}
  	if(ns==1) pv = sum(L[1,,TT])
	else pv = rowSums(L[,,TT])
  	lk = sum(yv*log(pv))
  	out = list(lk=lk,Phi=Phi,L=L,pv=pv)
}
