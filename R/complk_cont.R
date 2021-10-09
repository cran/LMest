complk_cont <- function(Y,piv,Pi,Mu,Si,k){

# Preliminaries
  	sY = dim(Y)
  	n = sY[1]
  	TT = sY[2]
  	if(length(sY)==2) r = 1 else r = sY[3]
  	if(r==1){
  		if(is.matrix(Y)) Y = array(Y,c(dim(Y),1))
#  		if(is.matrix(R)) R = array(R,c(dim(R),1))
  	}

# Compute log-likelihood
	Phi = array(1,c(n,k,TT)); L = array(0,c(n,k,TT))
	for(u in 1:k) Phi[,u,1] =  dmvnorm(matrix(Y[,1,],n,r),Mu[,u],Si)
	  	
  	L[,,1] = Phi[,,1]%*%diag(piv)
  	for(t in 2:TT){
		for(u in 1:k) Phi[,u,t] =  dmvnorm(matrix(Y[,t,],n,r),Mu[,u],Si)

   		L[,,t] = Phi[,,t]*(L[,,t-1]%*%Pi[,,t])
  	}
  	if(n==1) pv = sum(L[1,,TT])
	else pv = rowSums(L[,,TT])
  	lk = sum(log(pv))
  	out = list(lk=lk,Phi=Phi,L=L,pv=pv)
}
