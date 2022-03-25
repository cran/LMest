lk_comp_latent_cont <- function(Y,yv,Piv,PI,Mu,Si,k){

# der = TRUE for derivative
# Preliminaries
   	sY = dim(Y)
  	n = sY[1]
  	TT = sY[2]
  	if(length(sY)==2) r = 1 else r = sY[3]
  	if(r==1) if(is.matrix(Y)) Y = array(Y,c(dim(Y),1))
# Compute log-likelihood
	Phi = array(1,c(n,k,TT)); L = array(0,c(n,k,TT))
	
	for(u in 1:k) Phi[,u,1] = dmvnorm(matrix(Y[,1,],n,r),Mu[,u],Si)
#	Phi = pmax(Phi,10^-30)
	L[,,1] = Phi[,,1]*Piv
	for(t in 2:TT){
		for(u in 1:k) Phi[,u,t] = dmvnorm(matrix(Y[,t,],n,r),Mu[,u],Si)
#		Phi = pmax(Phi,10^-30)
   		for(i in 1:n)	L[i,,t] = L[i,,t-1]%*%PI[,,i,t]
   		L[,,t] = L[,,t]*Phi[,,t]
  	}
  	
  	if(n==1) pv = sum(L[,,TT])
  	else pv = rowSums(L[,,TT])
  	lk = sum(yv*log(pv))

# output
  	out = list(lk=lk,Phi=Phi,L=L,pv=pv)
}
