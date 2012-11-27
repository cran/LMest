lk_obs <-
function(th,Am,Bm,Cm,b,k,S,yv,T,r,mod){
	
# copute corresponding parameters
    ns = dim(S)[1]
    n = sum(yv)
	Psi = array(0,c(b+1,k,1))
	dim(Psi)
	for(u in 1:k){
		th1 = th[1:b]; th = th[-(1:b)]
		Psi[,u,1] = exp(Am%*%th1); Psi[,u,1] = Psi[,u,1]/sum(Psi[,u,1])				
	}
	th1 = th[1:(k-1)]; th = th[-(1:(k-1))]
	piv = exp(Bm%*%th1); piv = as.vector(piv/sum(piv))
	Pi = matrix(0,k,k)
	for(u in 1:k){
		th1 = th[1:(k-1)]; th = th[-(1:(k-1))]
		if(k==2) Pi[u,] = exp(Cm[,,u]*th1)
		else Pi[u,] = exp(Cm[,,u]%*%th1)
		Pi[u,] = Pi[u,]/sum(Pi[u,])
	}
	Pi = array(Pi,c(k,k,T))
	Pi[,,1] = 0
# compute log-likelihood
    out = complk(S,yv,piv,Pi,Psi,k)
  	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
  	sc = NULL
# ---- E-step ----
# Compute V and U
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
# ---- M-step ----
# Update Psi
   	Y = array(0,c(b+1,T,k,r))
   	for(j in 1:r) for(t in 1:T) for(jb in 0:b){
    	ind = which(S[,t,j]==jb)
   		li = length(ind)
    	if(li==1) Y[jb+1,t,,j] = V[ind,,t]
		if(li>1) Y[jb+1,t,,j] = colSums(V[ind,,t])
	}
 	Y1 = apply(Y,c(1,3,4),sum)
 	if(r==1) for(u in 1:k) sc = c(sc,t(Am)%*%(Y1[,u,1]-sum(Y1[,u,1])*Psi[,u,1]))
# Update piv and Pi
	sc = c(sc,t(Bm)%*%(colSums(V[,,1])-n*piv))
	if(mod==0) for(t in 2:T) Pi[,,t] = diag(1/rowSums(U[,,t]))%*%U[,,t]
	if(mod==1){
	   	Ut = apply(U[,,2:T],c(1,2),sum)
   	   	for(u in 1:k) sc = c(sc,t(Cm[,,u])%*%(Ut[u,]-sum(Ut[u,])*Pi[u,,2]))
	}
	if(mod==2){
#       	Ut1 = sum(U[,,2:mod],3);  #to correct
#	   	Ut2 = sum(U[,,mod+1:T],3);    #to correct
#       	Pi[,,2:mod] = repmat(diagv(1/sum(Ut1,2),Ut1),c(1,1,mod-1));       #to correct  
#       	Pi[,,mod+1:T] = repmat(diagv(1/sum(Ut2,2),Ut2),c(1,1,T-mod));     #to correct    
	}
# return
	out = list(lk=lk,sc=sc)	
}
