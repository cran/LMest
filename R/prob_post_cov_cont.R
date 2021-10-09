prob_post_cov_cont <- function(Y,yv,Mu,Si,Piv,PI,Phi,L,pv,der=FALSE,
                          dlPhi=NULL,dlPiv=NULL,dlPI=NULL,dlL=NULL,dlL2=NULL,dlpv=NULL){
	
# prelimiaries
   	sY = dim(Y)
  	n = sY[1]
  	TT = sY[2]
  	if(length(sY)==2) r = 1	else r = sY[3]
    k = ncol(Piv)
    
# use backward recursion for potesterior probabilities	
   	V = array(0,c(n,k,TT)); U = array(0,c(k,k,n,TT))
   	V1 = V
   	Pv = matrix(1/pv,n,k)
	Yv = matrix(yv,n,k); Yvp = Yv*Pv   		
   	# last step
	V1[,,TT] = L[,,TT]*Pv; V[,,TT] = Yv*V1[,,TT]
	AA = Yvp*Phi[,,TT]
	for(i in 1:n) U[,,i,TT] = outer(L[i,,TT-1],AA[i,])
	M = matrix(1,n,k)
	# backward
	if(TT>2) for(t in seq(TT-1,2,-1)){
		MP = Phi[,,t+1]*M
		for(i in 1:n) M[i,] = PI[,,i,t+1]%*%MP[i,]
		AA = Yvp*Phi[,,t]*M
		for(i in 1:n) U[,,i,t] = outer(L[i,,t-1],AA[i,])
   		V1[,,t] = L[,,t]*M*Pv; V[,,t] = Yv*V1[,,t]
    }
    U = U*PI
    MP = Phi[,,2]*M
    for(i in 1:n) M[i,] = PI[,,i,2]%*%MP[i,]
    V1[,,1] = L[,,1]*M*Pv; V[,,1] = Yv*V1[,,1]
       
# compute derivarives
# if the derivative is required
  	if(der){
  		ns = n #to check
  		nal = dim(dlPhi)[4]; nbe = dim(dlPiv)[3]; nga = dim(dlPI)[5]
  		npar = nal+nbe+nga
  		indal = 1:nal; indbe = nal+(1:nbe); indga = nal+nbe+(1:nga)
    	dlV = array(0,c(ns,k,TT,npar)); dlU = array(0,c(k,k,ns,TT,npar))
# last time
  		M = matrix(1,ns,k); dlM = array(0,c(ns,k,npar))
  		for(u in 1:k) dlV[,u,TT,] = dlL[,u,TT,] - dlpv
  		for(u in 1:k) for(ub in 1:k) dlU[ub,u,,TT,] = dlL2[,ub,u,TT,] - dlpv
# previous occasions
  		for(t in seq(TT-1,1,-1)){
  			M2 = array(0,c(ns,k,k)); dlM2 = array(0,c(ns,k,k,npar))
  			for(ub in 1:k) for(u in 1:k){
  				M2[,ub,u] = M[,u]*PI[ub,u,,t+1]*Phi[,u,t+1]
  				dlM2[,ub,u,] = dlM[,u,]
  				dlM2[,ub,u,indal] = dlM2[,ub,u,indal]+dlPhi[,u,t+1,]
  				dlM2[,ub,u,indga] = dlM2[,ub,u,indga]+dlPI[ub,u,,t+1,]
  			}
  			M = apply(M2,c(1,2),sum)
 		    dlM = array(0,c(ns,k,npar))
			for(ub in 1:k) for(u in 1:k) dlM[,ub,] = dlM[,ub,]+M2[,ub,u]/M[,ub]*dlM2[,ub,u,] 
			for(u in 1:k) dlV[,u,t,] = dlL[,u,t,] + dlM[,u,] - dlpv 			
 		    if(t>1) for(u in 1:k) for(ub in 1:k) dlU[ub,u,,t,] = dlL2[,ub,u,t,]+dlM[,u,] - dlpv
  		}
  	}else{
  		dlU = NULL; dlV=NULL
  	}
# final output
    out = list(U=U,V=V,dlU=dlU,dlV=dlV)
}
