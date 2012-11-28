est_lm_basic <-
function(S,yv,k,start=0,mod=0,tol=10^-6){

# Preliminaries
    Info = FALSE  # for the information matrix (to be implemented)
	n = sum(yv)
   	sS = dim(S)
  	ns = sS[1]
  	T = sS[2]
  	if(length(sS)==2) r = 1
  	else r = sS[3]
	Sv = matrix(S,ns*T,r)
  	b = max(S)
    Co = cbind(-diag(b),diag(b))
  	Ma = cbind(lower.tri(matrix(1,b,b), diag = TRUE),rep(0,b))
  	Ma = rbind(Ma,1-Ma)
  	th = NULL; sc = NULL
	J = NULL  		
  	if(Info){
  		A = cbind(-rep(1,b),diag(b))
  		Am = rbind(rep(0,b),diag(b))
  		B = cbind(-rep(1,k-1),diag(k-1))
  		Bm = rbind(rep(0,k-1),diag(k-1))
  		C = array(0,c(k-1,k,k))
  		Cm = array(0,c(k,k-1,k))
  		for(u in 1:k){
  			C[,,u] = rbind(cbind(diag(u-1),-rep(1,u-1),matrix(0,u-1,k-u)),
  						   cbind(matrix(0,k-u,u-1),-rep(1,k-u),diag(k-u)))
  			Cm[,,u] = rbind(cbind(diag(u-1),matrix(0,u-1,k-u)),
  							rep(0,k-1),
	 					    cbind(matrix(0,k-u,u-1),diag(k-u)))
  		}
  	}
# When there is just 1 latent class
  	if(k == 1){
	    piv = 1; Pi = 1
   	 	P = matrix(0,b+1,r)
	    for(t in 1:T){
	      	for(j in 1:r){
		        	for(y in 0:b){
		        		if(r==1) ind = which(S[,t]) else ind = which(S[,t,j]==y)
		  		        P[y+1,j] = P[y+1,j]+sum(yv[ind])
		        }
	   		}
	    }
	    Psi = P/(n*T)
	    pm = rep(1,ns)
	    for(t in 1:T) for(j in 1:r) pm = pm*Psi[S[,t,j]+1,j]
	    lk = sum(yv*log(pm))
	    np = r*b
	    aic = -2*lk+np*2
	    bic = -2*lk+np*log(n)
    	out = list(lk=lk,piv=piv,Pi=Pi,Psi=Psi,np=np,aic=aic,bic=bic,lkv=NULL,J=NULL,V=NULL,th=NULL,sc=NULL)
	    return(out)
  	}
# Starting values
	if(start == 0){
   		P = matrix(0,b+1,r)
        	for(t in 1:T){
        		for(j in 1:r){
           			for(y in 0:b){
	           			if(r==1) ind = which(S[,t]==y) else ind = which(S[,t,j]==y)
        					P[y+1,j] = P[y+1,j]+sum(yv[ind])
           			}
           		}
       		}
		E = Co%*%log(Ma%*%P)
  	   	Psi = array(0,c(b+1,k,r)); Eta = array(0,c(b,k,r))
        	grid = seq(-k,k,2*k/(k-1))
        	for(c in 1:k){
         		for(j in 1:r){
      			etac = E[,j]+grid[c]
          			Eta[,c,j] = etac
           			Psi[,c,j] = invglob(etac)
          		}
       	}
  		piv = rep(1,k)/k
    		Pi = matrix(1,k,k)+9*diag(k); Pi = diag(1/rowSums(Pi))%*%Pi;
    		Pi = array(Pi,c(k,k,T)); Pi[,,1] = 0
  	}
  	if(start==1){
	    Psi = array(runif((b+1)*k*r),c(b+1,k,r)) 
	    for(j in 1:r) for(c in 1:k) Psi[,c,j] = Psi[,c,j]/sum(Psi[,c,j])
	    Pi = array(runif(k^2*T),c(k,k,T))
	    for(t in 2:T) Pi[,,t] = diag(1/rowSums(Pi[,,t]))%*%Pi[,,t]
	    Pi[,,1] = 0
	    piv = runif(k); piv = piv/sum(piv)
	}
# Compute log-likelihood
    out = complk(S,yv,piv,Pi,Psi,k)
  	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
  	lk0 = sum(yv*log(yv/n)); dev = 2*(lk0-lk)
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
    cat("     mod    |      k      |    start    |     step    |     lk      |    lk-lko   | discrepancy |\n");
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
  	cat(sprintf("%11g",c(mod,k,start,0,lk)),"\n",sep=" | ")
  	it = 0; lko = lk-10^10; lkv = NULL
  	par = c(piv,as.vector(Pi),as.vector(Psi)); paro = par
# Iterate until convergence
	while(((abs(lk-lko)>tol | max(abs(par-paro))>tol) & it<10^5) || it<2){
		Psi0 = Psi; piv0 = piv; Pi0 = Pi
		it = it+1;
# ---- E-step ----
# Compute V and U
#time = proc.time()
   		V = array(0,c(ns,k,T)); U = array(0,c(k,k,T))
   		Yvp = matrix(yv/pv,ns,k)
  		M = matrix(1,ns,k);
   		V[,,T] = Yvp*L[,,T]
   		U[,,T] = (t(L[,,T-1])%*%(Yvp*Phi[,,T]))*Pi[,,T]
   		for(t in seq(T-1,2,-1)){
   			M = (Phi[,,t+1]*M)%*%t(Pi[,,t+1]);
      		V[,,t] = Yvp*L[,,t]*M
      		U[,,t] = (t(L[,,t-1])%*%(Yvp*Phi[,,t]*M))*Pi[,,t]
    	}
    	M = (Phi[,,2]*M)%*%t(Pi[,,2])
    	V[,,1] = Yvp*L[,,1]*M
#print(proc.time()-time)
# If required store parameters
# ---- M-step ----
# Update Psi
	Y1 = array(0,c(b+1,k,r))
	Vv = matrix(aperm(V,c(1,3,2)),ns*T,k)
	for(j in 1:r) for(jb in 0:b) {
		ind = which(Sv[,j]==jb)
		Y1[jb+1,,j] = colSums(Vv[ind,])				
	}
	for(j in 1:r) for(c in 1:k) Psi[,c,j] = Y1[,c,j]/sum(Y1[,c,j]) 
#print(proc.time()-time)
# Update piv and Pi
	piv = colSums(V[,,1])/n
	U = pmax(U,10^-300)
	if(mod==0) for(t in 2:T) Pi[,,t] = diag(1/rowSums(U[,,t]))%*%U[,,t]
	if(mod==1){
	    	Ut = apply(U[,,2:T],c(1,2),sum)
  	    	Pi[,,2:T] = array(diag(1/rowSums(Ut))%*%Ut,c(k,k,T-1))
    	}
    	if(mod>1){
	    	Ut1 = U[,,2:mod]
        	if(length(dim(Ut1))>2) Ut1 = apply(Ut1,c(1,2),sum)
        	Ut2 = U[,,(mod+1):T]
	    if(length(dim(Ut2))>2) Ut2 = apply(Ut2,c(1,2),sum)
    	   	Pi[,,2:mod] = array(diag(1/rowSums(Ut1,2))%*%Ut1,c(k,k,mod-1))         
        Pi[,,(mod+1):T] = array(diag(1/rowSums(Ut2,2))%*%Ut2,c(k,k,T-mod))         
	}
#print(proc.time()-time)
# Compute log-likelihood
    	paro = par; par = c(piv,as.vector(Pi),as.vector(Psi));
    	lko = lk;
    	out = complk(S,yv,piv,Pi,Psi,k)
    	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
    	if(it/10 == round(it/10)) cat(sprintf("%11g",c(mod,k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
    	lkv = c(lkv,lk)
#print(proc.time()-time)
	}
# Compute information matrix if required	
	if(Info){
		th = NULL
		for(u in 1:k) th = c(th,A%*%log(Psi[,u,1]))
		th = c(th,B%*%log(piv))
		for(u in 1:k) th = c(th,C[,,u]%*%log(Pi[u,,2]))
		
		lth = length(th)
		out = recursions(S,yv,Psi,piv,Pi,k,lth,Am,Bm,Cm,b)
		F1 = out$F1; F2 = out$F2; F1d = out$F1d; F2d = out$F2d

		sc = NULL
   		Y = array(0,c(b+1,T,k,r))
   		for(j in 1:r) for(t in 1:T) for(jb in 0:b){
      		ind = which(S[,t,j]==jb)
      		Y[jb+1,t,,j] = F1[,t,ind]%*%yv[ind]
   		}
   		Y1 = apply(Y,c(1,3,4),sum)
   		if(r==1){
   			for(u in 1:k) sc = c(sc,t(Am)%*%(Y1[,u,1]-sum(Y1[,u,1])*Psi[,u,1]))
   			J = sum(Y1[,1,1])*t(Am)%*%(diag(Psi[,1,1])-Psi[,1,1]%o%Psi[,1,1])%*%Am
   			for(u in 2:k) J = bdiag(J,sum(Y1[,u,1])*t(Am)%*%(diag(Psi[,u,1])-Psi[,u,1]%o%Psi[,u,1])%*%Am)
   		}
   	   	bv = F1[,1,]%*%yv
		sc = c(sc,t(Bm)%*%(bv-sum(bv)*piv))
		J = bdiag(J,n*t(Bm)%*%(diag(piv)-piv%o%piv)%*%Bm)
		Ut = 0
		for(i in 1:ns) for(t in 2:T) Ut = Ut+yv[i]*F2[,,t,i]
		for(u in 1:k){
   			sc = c(sc,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,2]))
   			J = bdiag(J,sum(Ut[u,])*t(Cm[,,u])%*%(diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2])%*%Cm[,,u])
   		}
   		J = as.matrix(J) 
   	   	
        Jd = NULL
        for(pa in 1:lth){
        	scj = NULL
	   		Y = array(0,c(b+1,T,k,r))
   			for(j in 1:r) for(t in 1:T) for(jb in 0:b){
   		   		ind = which(S[,t,j]==jb)
	   	   		Y[jb+1,t,,j] = F1d[,t,ind,pa]%*%yv[ind]
   			}
	   		Y1 = apply(Y,c(1,3,4),sum)
	   		if(r==1) for(u in 1:k) scj = c(scj,t(Am)%*%(Y1[,u,1]-sum(Y1[,u,1])*Psi[,u,1]))
	   	   	bv = F1d[,1,,pa]%*%yv
			scj = c(scj,t(Bm)%*%(bv-sum(bv)*piv))
			Ut = 0
			for(i in 1:ns) for(t in 2:T) Ut = Ut+yv[i]*F2d[,,t,i,pa]
			for(u in 1:k) scj = c(scj,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,2]))
			Jd = cbind(Jd,scj)
		}
		
		J = J-Jd
		V = ginv(J)
		for(u in 1:k){
			Om = diag(Psi[,u,1])-Psi[,u,1]%o%Psi[,u,1]
			if(u==1) M = Om%*%Am
			else M = bdiag(M,Om%*%Am)
		}
		Om = diag(piv)-tcrossprod(piv,piv)
		M = bdiag(M,Om%*%Bm)
		for(u in 1:k){
			Om = diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2]
			M = bdiag(M,Om%*%Cm[,,u])
		}
		M = as.matrix(M)
		V = M%*%V%*%t(M)

# to check derivatives
#		th0 = th-10^-5/2
#		out = lk_obs(th0,Am,Bm,Cm,b,k,S,yv,T,r,mod)
#		lk0 = out$lk; sc0 = out$sc
#		lth = length(th)
#		scn = rep(0,lth)
#		J = matrix(0,lth,lth)
#		for(j in 1:lth){
#			thj = th0; thj[j] = thj[j]+10^-5
#			out = lk_obs(thj,Am,Bm,Cm,b,k,S,yv,T,r,mod)
#			scn[j] = (out$lk-lk0)/10^-5
#			J[,j] = (out$sc-sc0)/10^-5
#		}
#		J = -(J+t(J))/2
#		print(cbind(sc,scn))
#		print(J)
#		print(J-J0)
	}
# Compute number of parameters  
	np = (k-1)+k*b*r
  	if(mod==0) np = np+(T-1)*k*(k-1)
  	if(mod==1) np = np+k*(k-1)
  	if(mod==2) np = np+2*k*(k-1)
  	aic = -2*lk+np*2
  	bic = -2*lk+np*log(n)
	cat(sprintf("%11g",c(mod,k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")	
	out = list(lk=lk,piv=piv,Pi=Pi,Psi=Psi,np=np,aic=aic,bic=bic,lkv=lkv,J=J,V=V,th=th,sc=sc)
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
	return(out)
}
