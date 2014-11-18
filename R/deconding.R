decoding <- function(est,Y,X1=NULL,X2=NULL,fort=TRUE){
	
	# Provide local decoding on the basis of the output of
	# est_lm_basic, est_lm_cov_latent, est_lm_cov_manifest (to be done, in this case the covariates are only in X1)
	# and est_lm_mixed
	
	# est = output from one of these function
	# Y = vector of responses or matrix of responses for which having local deconding
	#
	# Ul = matrix of local deconding corresponding to each row of Y
	# Ug = matrix of global deconding corresponding to each row of Y
	
# est_lm_basic	   
	if(class(est)=="LMbasic"){
		if(dim(est$Psi)[3]==1){
			if(is.vector(Y)) Y = t(Y)		
			n = nrow(Y); TT = ncol(Y)
		}else{
			if(is.matrix(Y)) Y = array(Y,c(1,dim(Y)))
			n = dim(Y)[1]; TT = dim(Y)[2]; r = dim(Y)[3]
		}
		piv = est$piv; Pi = est$Pi; Psi = est$Psi
		k = length(est$piv)
	    out = complk(Y,rep(1,n),piv,Pi,Psi,k)
	    Phi = out$Phi; L = out$L; pv = out$pv
   		V = array(0,c(n,k,TT))
   		Yvp = matrix(1/pv,n,k)
  		M = matrix(1,n,k)
   		V[,,TT] = Yvp*L[,,TT]
   		for(t in seq(TT-1,2,-1)){
   			M = (Phi[,,t+1]*M)%*%t(Pi[,,t+1])
      		V[,,t] = Yvp*L[,,t]*M
    	}
    	M = (Phi[,,2]*M)%*%t(Pi[,,2])
    	V[,,1] = Yvp*L[,,1]*M
# local deconding
   		Ul = matrix(0,n,TT)
		for(i in 1:n) for(t in 1:TT) Ul[i,t] = which.max(V[i,,t])
    	if(n==1) Ul = as.vector(Ul)
# global deconding (Viterbi)
		R = L; Ug = matrix(0,n,TT)
		for(i in 1:n) for(t in 2:TT) for(u in 1:k) R[i,u,t] = Phi[i,u,t]*max(R[i,,t-1]*Pi[,u,t])
		if(n==1) Ug[,TT] = which.max(R[,,TT])
		else Ug[,TT] = apply(R[,,TT],1,which.max)
		for(i in 1:n) for(t in seq(TT-1,1,-1)) Ug[i,t] = which.max(R[i,,t]*Pi[,Ug[i,t+1],t+1])
    	if(n==1) Ug = as.vector(Ug)
	}
	
# est_lm_cov_latent
	if(class(est)=="LMlatent"){
		param = est$param
		if(dim(est$Psi)[3]==1){
			if(is.vector(Y)) Y = t(Y)
			if(is.vector(X1)) X1 = t(X1)		
			if(is.matrix(X2)) X2 = array(X2,c(1,dim(X2)))		
			n = nrow(Y); TT = ncol(Y)
		}else{
			if(is.matrix(Y)) Y = array(Y,c(1,dim(Y)))
			if(is.vector(X1)) X1 = t(X1)		
			if(is.matrix(X2)) X2 = array(X2,c(1,dim(X2)))		
			n = dim(Y)[1]; TT = dim(Y)[2]; r = dim(Y)[3]
		}		
		k = ncol(est$Be)+1
		Psi = est$Psi
		nc1 = dim(X1)[2] # number of covariates on the initial probabilities		
		Xlab = 1:n
		if(k == 2){
			GBe = as.matrix(c(0,1))
		}else{
			GBe = diag(k); GBe = GBe[,-1]
		}
		XXdis = array(0,c(k,(k-1)*(nc1+1),n))
		for(i in 1:n){
			xdis = c(1,X1[i,])
			XXdis[,,i] = GBe%*%(diag(k-1)%x%t(xdis))
		}
		be = as.vector(est$Be)
		out = prob_multilogit(XXdis,be,Xlab,fort)
		Piv = out$P
		nc2 = dim(X2)[3] # number of covariates on the transition probabilities
		Z = NULL
		for(t in 1:(TT-1)) Z = rbind(Z,X2[,t,])
		Zlab = 1:(n*(TT-1)); Zndis = n*(TT-1)
		if(param=="multilogit"){
   	 	ZZdis = array(0,c(k,(k-1)*(nc2+1),Zndis,k))
		    for(h in 1:k){
			    if(k==2){
				    if(h == 1) GGa = as.matrix(c(0,1)) else GGa = as.matrix(c(1,0))
			    }else{
				    GGa = diag(k); GGa = GGa[,-h]
			    }  		
			    for(i in 1:Zndis){
				    zdis = c(1,Z[i,])
				    ZZdis[,,i,h] = GGa%*%(diag(k-1)%x%t(zdis))
			    }
		    }
		 }else if(param=="difflogit"){
   	     Zlab = (((Zlab-1)*k)%x%rep(1,k))+rep(1,n*(TT-1))%x%(1:k)
   	     ZZdis = array(0,c(k,k*(k-1)+(k-1)*nc2,Zndis*k))
   	     j = 0
			for(i in 1:Zndis){
   	         for(h in 1:k){
   	             j = j+1
   	             if(k==2){
				      if(h == 1) GGa = as.matrix(c(0,1)) else GGa = as.matrix(c(1,0))
			        }else{
				        GGa = diag(k); GGa = GGa[,-h]
			        }  		
				    u = matrix(0,1,k); u[1,h] = 1
				    U = diag(k); U[,h] = U[,h]-1
				    U = U[,-1]
			        ZZdis[,,j] = cbind(u%x%GGa,U%x%t(Z[i,]))            
   	         	}
		    }
		}
        if(param=="multilogit"){
    		Ga = matrix(est$Ga,(nc2+1)*(k-1),k)
	    	PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,n,TT))
		    for(h in 1:k){
			    out = prob_multilogit(ZZdis[,,,h],Ga[,h],Zlab,fort)
			    PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,n,TT-1))
		    }
		}else if(param=="difflogit"){
            Ga = c(as.vector(t(est$Ga[[1]])),as.vector(est$Ga[[2]]))
            PI = array(0,c(k,k,n,TT))
            out = prob_multilogit(ZZdis,Ga,Zlab,fort)
   		    PIdis = out$Pdis; PI[,,,2:TT] = array(as.vector(t(out$P)),c(1,k,n,TT-1))
		}
		out = lk_comp_latent(Y,rep(1,n),Piv,PI,Psi,k,fort=fort)
		Phi = out$Phi; L = out$L; pv = out$pv
		out = prob_post_cov(Y,rep(1,n),Psi,Piv,PI,Phi,L,pv,fort=fort)
		V = out$V
# local deconding
   		Ul = matrix(0,n,TT)
		for(i in 1:n) for(t in 1:TT) Ul[i,t] = which.max(V[i,,t])
    	if(n==1) Ul = as.vector(Ul)
# global deconding (Viterbi)
		R = L; Ug = matrix(0,n,TT)
		for(i in 1:n) for(t in 2:TT) for(u in 1:k) R[i,u,t] = Phi[i,u,t]*max(R[i,,t-1]*PI[,u,i,t])
		if(n==1) Ug[,TT] = which.max(R[,,TT])
		else Ug[,TT] = apply(R[,,TT],1,which.max)
		for(i in 1:n) for(t in seq(TT-1,1,-1)) Ug[i,t] = which.max(R[i,,t]*PI[,Ug[i,t+1],i,t+1])
    	if(n==1) Ug = as.vector(Ug)		
	}

# est_lm_mixed
	if(class(est)=="LMmixed"){
		if(dim(est$Psi)[3]==1){
			if(is.vector(Y)) Y = t(Y)		
			n = nrow(Y); TT = ncol(Y)
		}else{
			if(is.matrix(Y)) Y = array(Y,c(1,dim(Y)))
			n = dim(Y)[1]; TT = dim(Y)[2]; r = dim(Y)[3]
		}
		yv = rep(1,n)
		la = est$la; Piv = est$Piv; Pi = est$Pi; Psi = est$Psi
		k1 = length(la); k2 = nrow(Piv)		
		Fc1 = matrix(0,n,k1); Fc2 = array(0,c(n,k1,k2)); Fc3 = array(0,c(n,k1,k2,k2))
		PP1 = array(0,c(n,k1,k2,TT))
		Phi = array(1,c(n,k2,TT))
		for(t in 1:TT){
  			if(r==1) Phi[,,t] = Phi[,,t]*Psi[Y[,t]+1,,1] 
  			else for(j in 1:r) Phi[,,t] = Phi[,,t]*Psi[Y[,t,j]+1,,j]
		}
		for(i in 1:n) for(u in 1:k1){
	   		 o = .Fortran("BWforback", TT, k2, Phi[i,,], Piv[,u], Pi[,,u], lk=0, Pp1=matrix(0,k2,TT), 
	        	         Pp2=array(0,c(k2,k2,TT)))
			Fc1[i,u] = exp(o$lk); Fc2[i,u,] = o$Pp1[,1]; Fc3[i,u,,] = apply(o$Pp2,c(1,2),sum)
  			PP1[i,u,,] = o$Pp1
  		}
		Fj1 = Fc1%*%diag(la)
		fm = rowSums(Fj1)
		fm = pmax(fm,10^-300)
		W = (Fj1/matrix(fm,n,k1))*yv
			U1 = apply(W,1,which.max)
		PV = array(W,c(n,k1,k2,TT))*PP1
# local deconding
   		Ul = matrix(0,n,TT)
		for(i in 1:n) for(t in 1:TT) Ul[i,t] = which.max(PV[i,U1[i],,t])
    	if(n==1) Ul = as.vector(Ul)
# global deconding (Viterbi)
		R = array(0,c(n,k2,TT))
		for(i in 1:n) for(v in 1:k2) R[i,v,1] = Phi[i,v,1]*Piv[v,U1[i]]		 
		for(i in 1:n) for(t in 2:TT) for(v in 1:k2) R[i,v,t] = Phi[i,v,t]*max(R[i,,t-1]*Pi[,v,U1[i]])
		Ug = matrix(0,n,TT)
		if(n==1) Ug[,TT] = which.max(R[,,TT])
		else Ug[,TT] = apply(R[,,TT],1,which.max)
		for(i in 1:n) for(t in seq(TT-1,1,-1)) Ug[i,t] = which.max(R[i,,t]*Pi[,Ug[i,t+1],U1[i]])
    	if(n==1) Ug = as.vector(Ug)
	}
	
#est_lm_manifest
if(class(est)=="LMmanifest"){
	if(is.vector(Y)) Y = t(Y)		
	n = nrow(Y); TT = ncol(Y)
	k = dim(est$PI)[2]
	V = est$PRED0
	Phi = est$Phi
	piv = est$la
	PI = est$PI
	# local deconding
	Ul = matrix(0,n,TT)
	for(i in 1:n) for(t in 1:TT) Ul[i,t] = which.max(V[i,,t])
 	if(n==1) Ul = as.vector(Ul)
	# global deconding (Viterbi)
	R = array(0,c(n,k,TT))
	for(i in 1:n) for(v in 1:k) R[i,v,1] = Phi[i,v,1]*piv[v]	
	Ug = matrix(0,n,TT)
	for(i in 1:n) for(t in 2:TT) for(u in 1:k) R[i,u,t] = Phi[i,u,t]*max(R[i,,t-1]*PI[,u])
	if(n==1) Ug[,TT] = which.max(R[,,TT])
	else Ug[,TT] = apply(R[,,TT],1,which.max)
	for(i in 1:n) for(t in seq(TT-1,1,-1)) Ug[i,t] = which.max(R[i,,t]*PI[,Ug[i,t+1]])
	if(n==1) Ug = as.vector(Ug)	
}
# output
	out = list(Ul=Ul,Ug=Ug)
	
}