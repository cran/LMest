lk_comp_latent <- function(S,yv,Piv,PI,Psi,k,fort=TRUE,der=FALSE,dlPsi=NULL,dlPiv=NULL,dlPI=NULL){

# der = TRUE for derivative

# Preliminaries
   	sS = dim(S)
  	ns = sS[1]
  	TT = sS[2]
  	l = max(S)+1
  	if(length(sS)==2) r = 1 else r = sS[3]
# Compute log-likelihood
	Phi = array(1,c(ns,k,TT)); L = array(0,c(ns,k,TT))
	if(fort){
        if(r==1){
	       for(t in 1:TT) Phi[,,t] = Phi[,,t]*Psi[S[,t]+1,,1]
	    }else{
            o = .Fortran("for_mult",as.integer(TT),as.integer(r),as.integer(k),as.integer(ns),as.integer(l),as.integer(S),Psi,Piv,PI,Phi=Phi,LL=array(0,c(ns,k,TT)))
            Phi = o$Phi; L = o$LL
	    }
	}else{
		if(r==1){
	        for(t in 1:TT) Phi[,,t] = Phi[,,t]*Psi[S[,t]+1,,1]
	    }else{
	        for(t in 1:TT) for(j in 1:r) Phi[,,t] = Phi[,,t]*Psi[S[,t,j]+1,,j]
	    }
	}  	
    if(r==1 || fort==F){
      	L[,,1] = Phi[,,1]*Piv
		for(t in 2:TT){	
  			for(i in 1:ns)	L[i,,t] = L[i,,t-1]%*%PI[,,i,t]
  			L[,,t] = L[,,t]*Phi[,,t]
  		}
  	}
  	if(ns==1) pv = sum(L[,,TT])
  	else pv = rowSums(L[,,TT])
  	lk = sum(yv*log(pv))
# if the derivative is required
  	if(der){
  		nal = dim(dlPsi)[4]; nbe = dim(dlPiv)[3]; nga = dim(dlPI)[5]
  		npar = nal+nbe+nga
  		indal = 1:nal; indbe = nal+(1:nbe); indga = nal+nbe+(1:nga)
        dlPhi = array(0,c(ns,k,TT,nal)); dlL = array(0,c(ns,k,TT,nal+nbe+nga))
		if(r==1){
	        for(t in 1:TT) dlPhi[,,t,] = dlPsi[S[,t]+1,,1,]
	    }else{
	        for(t in 1:TT) for(j in 1:r) dlPhi[,,t,] = dlPhi[,,t,]+dlPsi[S[,t,j]+1,,j,]
	    }
	    dlL[,,1,indal] = dlPhi[,,1,]; dlL[,,1,indbe] = dlPiv
		dlL2 = array(0,c(ns,k,k,t,npar))
		for(t in 2:TT){
			L2 = array(0,c(ns,k,k))
			for(ub in 1:k) for(u in 1:k){
   				L2[,ub,u] = L[,ub,t-1]*PI[ub,u,,t]*Phi[,u,t]
   				dlL2[,ub,u,t,] = dlL[,ub,t-1,]
   				dlL2[,ub,u,t,indal] = dlL2[,ub,u,t,indal]+dlPhi[,u,t,]
   				dlL2[,ub,u,t,indga] = dlL2[,ub,u,t,indga]+dlPI[ub,u,,t,]
  			}
            for(ub in 1:k) for(u in 1:k) dlL[,u,t,] = dlL[,u,t,]+L2[,ub,u]/L[,u,t]*dlL2[,ub,u,t,]
        }
        dlpv = matrix(0,ns,npar)
        for(u in 1:k) dlpv = dlpv+L[,u,TT]/pv*dlL[,u,TT,]
        dlk = colSums(dlpv)
  	}else{
  		dlk=NULL; dlL=NULL; dlL2=NULL; dlpv = NULL; dlPhi=NULL
  	}
# output
  	out = list(lk=lk,Phi=Phi,L=L,pv=pv,dlk=dlk,dlPhi=dlPhi,dlL=dlL,dlL2=dlL2,dlpv=dlpv)
}
