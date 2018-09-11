est_lm_basic_cont <-
function(Y,k,start=0,mod=0,tol=10^-8,maxit=1000,piv=NULL,Pi=NULL,Mu=NULL,Si=NULL){

# Preliminaries
    check_der = FALSE  # to check derivatives
   	sY = dim(Y)
  	n = sY[1]
  	TT = sY[2]
  	 
  	if(is.data.frame(Y)){
  		warning("Data frame not allowed for Y")
  	}
  	
  	if(length(sY)==2){
  		r = 1
  		if(is.matrix(Y)) Y = array(Y,c(dim(Y),1))
  	}else r = sY[3]
  	
  	# miss = any(is.na(S))
	# if(miss){
         # cat("Missing data in the dataset, treated as missing at random\n")
         # R = 1 * (!is.na(S))
         # S[is.na(S)] = 0
	# }else{
		# R = NULL
	# }
  	Yv = matrix(Y,n*TT,r)
  	
# When there is just 1 latent class
  	if(k == 1){ 
	    piv = 1; Pi = 1
   	 	Mu = colMeans(Yv)
   	 	Si = cov(Yv)
   	 	lk = sum(dmvnorm(Yv,Mu,Si,log=TRUE))
		np = k*r+r*(r+1)/2
	    aic = -2*lk+np*2
	    bic = -2*lk+np*log(n)
    		out =     		list(lk=lk,piv=piv,Pi=Pi,Mu=Mu,Si=Si,np=np,aic=aic,bic=bic,lkv=NULL,V=NULL,call=match.call())
		class(out)="LMbasiccont" 
	    return(out)
  	}
  	
# Starting values	
	if(start == 0){
		mu = colMeans(Yv)
		Si = cov(Yv); std = sqrt(diag(Si))
		qt = qnorm((1:k)/(k+1))
		Mu = matrix(0,r,k)
		for(u in 1:k) Mu[,u] = qt[u]*std+mu
   		
  		piv = rep(1,k)/k
		Pi = matrix(1,k,k)+9*diag(k); Pi = diag(1/rowSums(Pi))%*%Pi;
		Pi = array(Pi,c(k,k,TT)); Pi[,,1] = 0
  	}
  	if(start==1){
  		Mu = matrix(0,r,k)
		mu = colMeans(Yv)
		Si = cov(Yv)
		for(u in 1:k) Mu[,u] = rmvnorm(1,mu,Si)
	    Pi = array(runif(k^2*TT),c(k,k,TT))
	    for(t in 2:TT) Pi[,,t] = diag(1/rowSums(Pi[,,t]))%*%Pi[,,t]
	    Pi[,,1] = 0
	    piv = runif(k); piv = piv/sum(piv)
	}
	if(start==2){
		if(is.null(piv)) stop("initial value of the initial probabilities (piv) must be given in input")
		if(is.null(Pi)) stop("initial value of the transition probabilities (Pi) must be given in input")
		if(is.null(Mu)) stop("initial value of the conditional means of the response variables (Mu) must be given in input")
		if(is.null(Si)) stop("initial value of the var-cov matrix common to all states (Si) must be given in input")
		piv = piv
		Pi = Pi
		Mu = Mu 
		Si = Si
	}
	
# Compute log-likelihood
    out = complk_cont(Y,piv,Pi,Mu,Si,k)
  	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
    cat("     mod    |      k      |    start    |     step    |     lk      |    lk-lko   | discrepancy |\n");
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
  	cat(sprintf("%11g",c(mod,k,start,0,lk)),"\n",sep=" | ")
  	it = 0; lko = lk-10^10; lkv = NULL
  	par = c(piv,as.vector(Pi),as.vector(Mu),as.vector(Si))
  	if(any(is.na(par))) par = par[-which(is.na(par))]
  	paro = par
# Iterate until convergence
	while((lk-lko)/abs(lk)>tol & it<maxit){
		Mu0 = Mu; Si0 = Si; piv0 = piv; Pi0 = Pi
		it = it+1;
# ---- E-step ----
# Compute V and U
#time = proc.time()
   		V = array(0,c(n,k,TT)); U = array(0,c(k,k,TT))
   		Yvp = matrix(1/pv,n,k)
  		M = matrix(1,n,k)
   		V[,,TT] = Yvp*L[,,TT]
   		U[,,TT] = (t(L[,,TT-1])%*%(Yvp*Phi[,,TT]))*Pi[,,TT]
   		if(TT>2){
   			for(t in seq(TT-1,2,-1)){
   				M = (Phi[,,t+1]*M)%*%t(Pi[,,t+1]);
      			V[,,t] = Yvp*L[,,t]*M
      			U[,,t] = (t(L[,,t-1])%*%(Yvp*Phi[,,t]*M))*Pi[,,t]
			}
		}		
		M = (Phi[,,2]*M)%*%t(Pi[,,2])
		V[,,1] = Yvp*L[,,1]*M

# If required store parameters
# ---- M-step ----
# Update Mu 
	Vv = matrix(aperm(V,c(1,3,2)),n*TT,k)
 	for(u in 1:k) Mu[,u] = (t(Yv)%*%Vv[,u])/sum(Vv[,u])
# Update Si 
     Si = matrix(0,r,r)      		
    for(u in 1:k) Si= Si+ t(Yv-rep(1,n*TT)%*%t(Mu[,u]))%*%diag(Vv[,u])%*%as.matrix(Yv-rep(1,n*TT)%*%t(Mu[,u]))        
	Si = Si/(n*TT)

#print(proc.time()-time)
# Update piv and Pi
	piv = colSums(V[,,1])/n
	U = pmax(U,10^-300)
	if(mod==0) for(t in 2:TT) Pi[,,t] = diag(1/rowSums(U[,,t]))%*%U[,,t]
	if(mod==1){
	    	Ut = apply(U[,,2:TT],c(1,2),sum)
  	    	Pi[,,2:TT] = array(diag(1/rowSums(Ut))%*%Ut,c(k,k,TT-1))
    	}
    	if(mod>1){
	    	Ut1 = U[,,2:mod]
        	if(length(dim(Ut1))>2) Ut1 = apply(Ut1,c(1,2),sum)
        	Ut2 = U[,,(mod+1):TT]
	    if(length(dim(Ut2))>2) Ut2 = apply(Ut2,c(1,2),sum)
    	   	Pi[,,2:mod] = array(diag(1/rowSums(Ut1,2))%*%Ut1,c(k,k,mod-1))         
        Pi[,,(mod+1):TT] = array(diag(1/rowSums(Ut2,2))%*%Ut2,c(k,k,TT-mod))         
	}
#print(proc.time()-time)
# Compute log-likelihood
    	paro = par; par = c(piv,as.vector(Pi),as.vector(Mu),as.vector(Si))
    if(any(is.na(par))) par = par[-which(is.na(par))]
    	lko = lk
    	out = complk_cont(Y,piv,Pi,Mu,Si,k)
 
    	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
    	if(it/10 == round(it/10)) cat(sprintf("%11g",c(mod,k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
    	lkv = c(lkv,lk)
#print(proc.time()-time)
	}

# Compute number of parameters  
	np = (k-1)+k*r+r*(r+1)/2
  	if(mod==0) np = np+(TT-1)*k*(k-1)
  	if(mod==1) np = np+k*(k-1)
  	if(mod>1) np = np+2*k*(k-1)
  	aic = -2*lk+np*2
  	bic = -2*lk+np*log(n)
	cat(sprintf("%11g",c(mod,k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")	
	# adjust output
#	if(any(yv!=1)) V = V/yv
	
	lk = as.vector(lk)
	dimnames(Pi)=list(state=1:k,state=1:k,time=1:TT)
#	dimnames(Mu) = list(dimnames(Y)[[3]],state=1:k)
#	dimnames(Si) = list(dimnames(Y)[[3]],dimnames(Y)[[3]])
	if(r==1) dimnames(Mu) = list(item=1,state=1:k) else dimnames(Mu)=list(item=1:r,state=1:k)
	dimnames(Si)=list(item=1:r,item=1:r)

	out = list(lk=lk,piv=piv,Pi=Pi,Mu=Mu,Si=Si,np=np,aic=aic,bic=bic,lkv=lkv,V=V,call=match.call())
	
    cat("------------|-------------|-------------|-------------|-------------|-------------|-------------|\n");
    class(out)="LMbasiccont"
	return(out)
}
