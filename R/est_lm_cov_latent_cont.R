est_lm_cov_latent_cont <-
function(Y,X1=NULL,X2=NULL,yv = rep(1,nrow(Y)),k,start=0,tol=10^-8,maxit=1000,
         param="multilogit",Mu=NULL,Si=NULL,Be=NULL,Ga=NULL,output=FALSE){

# Fit the LM model for continuous outcomes with individual covariates in the distribution of the latent process
#
# INPUT:
# Y = array of available continuous outcome (n x TT x r)
# X1 = matrix of covariates affecting the initial probabilities
# X2 = array of covariates affecting the transition probabilities
# yv = vector of frequencies
# k = number of latent states
# start = initialization (0 = deterministic, 1 = random, 2 = initial values in input)
# maxit = maximum number of iterations
# param = type of parametrization for the transition probabilities:
#         multilogit = standard multinomial logit for every row of the transition matrix
#         difflogit  = multinomial logit based on the difference between two sets of parameters
# Mu = conditional means of the response variables (if start=2)
# Si = var-cov matrix common to all states (if start=2)
# Be = parameters on the initial probabilities (if start=2)
# Ga = parameters on the transition probabilities (if start=2)
# output = to return additional output

# Preliminaries
    check_der = FALSE # to check score and info
   	sY = dim(Y)
  	n = sY[1]
  	TT = sY[2]
  	if(length(sY)==2) r = 1
  	else r = sY[3]
  	if(is.data.frame(Y)) warning("Data frame not allowed for Y")
   	if(!is.null(X1)) if(any(is.na(X1))) stop("missing data not allowed in X1")
   	if(!is.null(X2)) if(any(is.na(X2))) stop("missing data not allowed in X2")
  	
	Yv = matrix(Y,n*TT,r)
	  
	## Check and inpute for missing data
	
	miss = any(is.na(Yv))
	if(miss)
	{
		Yv = cbind(1,Yv)
		pYv = prelim.mix(Yv,1)
		thhat = em.mix(prelim.mix(Yv,1))
		rngseed(1)
		Yv = imp.mix(pYv, da.mix(pYv,thhat,steps=100), Yv)[,-1]
		Y = array(Yv,c(n,TT,r))
		cat("Missing data in the dataset. imp.mix function (mix package) used for imputation.\n")
	}
			
# Covariate structure and related matrices: initial probabilities
	if(k == 2){
		GBe = as.matrix(c(0,1))
	}else{
		GBe = diag(k); GBe = GBe[,-1]
	}
	if(is.null(X1)){
		nc1=0
		Xlab = rep(1,n)
		nameBe = NULL
	}else{
		if(is.vector(X1)) X1 = matrix(X1,n,1)
		nc1 = dim(X1)[2] # number of covariates on the initial probabilities
		if(n!= dim(X1)[1]) stop("dimension mismatch between S and X1") 
		nameBe = colnames(X1)	
		out = aggr_data(X1)
		Xdis = out$data_dis
		if(nc1==1) Xdis = matrix(Xdis,length(Xdis),1)
		Xlab = out$label
	}	
	Xndis = max(Xlab)
	XXdis = array(0,c(k,(k-1)*(nc1+1),Xndis))
	for(i in 1:Xndis){
		if(nc1==0) xdis = 1 else xdis = c(1,Xdis[i,])
		XXdis[,,i] = GBe%*%(diag(k-1)%x%t(xdis))
	}
		

# for the transition probabilities
	if(is.null(X2)){
		if(param=="difflogit"){
			warning("with X2=NULL parametrization difflogit not considered")	
			param="multilogit"
		}
		nc2 = 0
		Zlab = rep(1,n*(TT-1))
		nameGa = NULL
		Zndis = max(Zlab)
	}else{
		if(TT==2) X2 = array(X2,c(n,1,dim(X2)[2]))
		if(is.matrix(X2)) X2 = array(X2,c(n,TT-1,1))
    	nc2 = dim(X2)[3] # number of covariates on the transition probabilities
    	if(n!= dim(X2)[1]) stop("dimension mismatch between S and X2")
    	nameGa = colnames(aperm(X2,c(1,3,2)))
		Z = NULL
		for(t in 1:(TT-1)) Z = rbind(Z,X2[,t,])
		if(nc2==1) Z = as.vector(X2)
		out = aggr_data(Z); Zdis = out$data_dis; Zlab = out$label; Zndis = max(Zlab)
		if(nc2==1) Zdis=matrix(Zdis,length(Zdis),1)
	}	
	if(param=="multilogit"){
    		ZZdis = array(0,c(k,(k-1)*(nc2+1),Zndis,k))
	    	for(h in 1:k){
			if(k==2){
				if(h == 1) GGa = as.matrix(c(0,1)) else GGa = as.matrix(c(1,0))
		    	}else{
			    	GGa = diag(k); GGa = GGa[,-h]
		    	}	  		
		    	for(i in 1:Zndis){
			    	if(nc2==0) zdis = 1 else zdis = c(1,Zdis[i,])
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
		        	ZZdis[,,j] = cbind(u%x%GGa,U%x%t(Zdis[i,]))            
            	}
	    	}
    	}
	
# When there is just 1 latent class
  	if(k == 1){
		Piv = rep(1,n); Pi = 1
		yvv = rep(yv,TT)
		Mu = colSums(yvv*Yv)/sum(yvv)
		Di = Yv-rep(1,n*TT)%o%Mu
		Si = t(Di)%*%(yvv*Di)/sum(yvv)
		lk = sum(yvv*dmvnorm(Yv,Mu,Si,log=TRUE))
		np = k*r+r*(r+1)/2
	    aic = -2*lk+np*2
	    bic = -2*lk+np*log(n)
		out = list(lk=lk,Piv=Piv,Pi=Pi,Mu=Mu,Si=Si,np=np,aic=aic,bic=bic,lkv=NULL,V=NULL,call=match.call())
		class(out)="LMlatentcont"
	    return(out)
  	}
# Starting values: deterministic initialization
	if(start == 0){
		yvv = rep(yv,TT)
		mu = colSums(yvv*Yv)/sum(yvv)
		Di = Yv-rep(1,n*TT)%o%mu
		Si = t(Di)%*%(yvv*Di)/sum(yvv)
		std = sqrt(diag(Si))
		qt = qnorm((1:k)/(k+1))
		Mu = matrix(0,r,k)
		for(u in 1:k) Mu[,u] = qt[u]*std+mu
# parameters on initial probabilities
       	be = array(0,(nc1+1)*(k-1))
       	out = prob_multilogit(XXdis,be,Xlab)
       	Piv = out$P; Pivdis = out$Pdis
# parameters on transition probabilities
        if(param=="multilogit"){
            Ga = matrix(0,(nc2+1)*(k-1),k)
            Ga[1+(0:(k-2))*(nc2+1),] = -log(10)
			PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,n,TT))
			for(h in 1:k){		
				tmp = ZZdis[,,,h]
				if(nc2==0) tmp = array(tmp,c(k,(k-1),Zndis))
				out = prob_multilogit(tmp,Ga[,h],Zlab)
			    PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,n,TT-1))
		    }		  
		}else if(param=="difflogit"){
             Ga = matrix(0,k*(k-1)+(k-1)*nc2)
             Ga[1:((h-1)*k)] = -log(10)
             PI = array(0,c(k,k,n,TT))
             out = prob_multilogit(ZZdis,Ga,Zlab)
   		     PIdis = out$Pdis; PI[,,,2:TT] = array(as.vector(t(out$P)),c(k,k,n,TT-1))
   		     PI = aperm(PI,c(2,1,3,4))
		}
  	}
  	
# random initialization
  	if(start==1){
		Mu = matrix(0,r,k)
		mu = colMeans(Yv)
		Si = cov(Yv)
		for(u in 1:k) Mu[,u] = rmvnorm(1,mu,Si)
# parameters on initial probabilities
       	be = c(rnorm(1),rep(0,nc1))
       	if(k>2) for(h in 2:(k-1)) be = c(be,rnorm(1),rep(0,nc1))
       	out = prob_multilogit(XXdis,be,Xlab)
       	Piv = out$P; Pivdis = out$Pdis
# parameters on transition probabilities
        if(param=="multilogit"){
    	#	Ga = matrix(-abs(rnorm((nc2+1)*(k-1),k)),(nc2+1)*(k-1),k)/2
    		Ga = matrix(0,(nc2+1)*(k-1),k)
    		Ga[1+(0:(k-2))*(nc2+1),] = -abs(rnorm((k-1)))   
	    	PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,n,TT))
		    for(h in 1:k){
		    		tmp = ZZdis[,,,h]
				if(nc2==0) tmp = array(tmp,c(k,(k-1),Zndis))
			    out = prob_multilogit(tmp,Ga[,h],Zlab)
			    PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,n,TT-1))
		   }
    	}else if(param=="difflogit"){
            Ga = c(-abs(rnorm(k*(k-1))),rep(0,(k-1)*nc2))
            PI = array(0,c(k,k,n,TT))
            out = prob_multilogit(ZZdis,Ga,Zlab)
   		    PIdis = out$Pdis; PI[,,,2:TT] = array(as.vector(t(out$P)),c(k,k,n,TT-1))
   		    PI = aperm(PI,c(2,1,3,4))
		}
	}
# initialization as input
    if(start==2){
# parameters on initial probabilities
       	be = as.vector(Be)
       	out = prob_multilogit(XXdis,be,Xlab)
       	Piv = out$P; Pivdis = out$Pdis
# parameters on transition probabilities
        if(param=="multilogit"){
        	if(is.list(Ga)) stop("invalid mode (list) for Ga")
    		Ga = matrix(Ga,(nc2+1)*(k-1),k)
	    	PIdis = array(0,c(Zndis,k,k)); PI = array(0,c(k,k,n,TT))
		    for(h in 1:k){
			    out = prob_multilogit(ZZdis[,,,h],Ga[,h],Zlab)
			    PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,n,TT-1))
		    }
		}else if(param=="difflogit"){
			if(is.list(Ga)) Ga = c(as.vector(t(Ga[[1]])),as.vector(Ga[[2]]))
            if(length(Ga)!=k*(k-1)+(k-1)*nc2) stop("invalid dimensions for Ga")
            PI = array(0,c(k,k,n,TT))
            out = prob_multilogit(ZZdis,Ga,Zlab)
   		    PIdis = out$Pdis; PI[,,,2:TT] = array(as.vector(t(out$P)),c(k,k,n,TT-1))
   		    PI = aperm(PI,c(2,1,3,4))
		}
    }

###### standard EM #####
   	out = lk_comp_latent_cont(Y,yv,Piv,PI,Mu,Si,k) 
   	lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
  # 	if(is.nan(lk)) browser()
	it = 0; lko = lk-10^10; lkv = NULL
	par = c(as.vector(Piv),as.vector(PI),as.vector(Mu),as.vector(Si))
	if(any(is.na(par))) par = par[-which(is.na(par))]
	paro = par
# Iterate until convergence
# display output	
	cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
    cat("      k     |    start    |     step    |     lk      |    lk-lko   | discrepancy |\n");
    cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
  	cat(sprintf("%11g",c(k,start,0,lk)),"\n",sep=" | ")
   	#cat("",sprintf("%11g",c(0,lk)),"\n",sep=" | ")
	while((lk-lko)/abs(lk)>tol & it<maxit){
		Mu0 = Mu; Si0 = Si; Piv0 = Piv; PI0 = PI
		it = it+1
# ---- E-step ----
# Compute V and U
       	out = prob_post_cov_cont(Y,yv,Mu,Si,Piv,PI,Phi,L,pv)
       	U = out$U; V = out$V  
# If required store parameters
# ---- M-step ----
# Update Mu
		Vv = matrix(aperm(V,c(1,3,2)),n*TT,k)
		for(u in 1:k) Mu[,u] = (t(Yv)%*%Vv[,u])/sum(Vv[,u])
# Update Si 
		Si = matrix(0,r,r)
		for(u in 1:k) Si= Si+ t(Yv-rep(1,n*TT)%o%Mu[,u])%*%diag(Vv[,u])%*%
		                      as.matrix(Yv-rep(1,n*TT)%o%Mu[,u])
		Si = Si/(sum(yv)*TT)
# Update piv
		out = est_multilogit(V[,,1],XXdis,Xlab,be,Pivdis)
		be = out$be; Pivdis = out$Pdi; Piv = out$P
# Update Pi
		if(param=="multilogit"){
	    	for(h in 1:k){
		    	UU = NULL
		    	for(t in 2:TT) UU = rbind(UU,t(U[h,,,t]))
		    	tmp = ZZdis[,,,h]
				if(nc2==0) tmp = array(tmp,c(k,(k-1),Zndis))
				tmp2 = PIdis[,,h]
				if(Zndis==1) tmp2 = matrix(tmp2,1,k)
	    		out = est_multilogit(UU,tmp,Zlab,Ga[,h],tmp2)
		    	PIdis[,,h] = out$Pdis; PI[h,,,2:TT] = array(as.vector(t(out$P)),c(1,k,n,TT-1)); Ga[,h] = out$be
	   	 	}
		}else if(param=="difflogit"){
		    Tmp = aperm(U[,,,2:TT],c(1,3,4,2))
		    Tmp = matrix(Tmp,n*k*(TT-1),k)
           	out = est_multilogit(Tmp,ZZdis,Zlab,Ga,PIdis)
	    		PIdis = out$Pdis; Ga = out$be
	    		Tmp = array(out$P,c(k,n,TT-1,k))
	    		PI[,,,2:TT] = aperm(Tmp,c(1,4,2,3)) 
        }
# Compute log-likelihood
   		paro = par; par = c(as.vector(Piv),as.vector(PI),as.vector(Mu),as.vector(Si))
   		if(any(is.na(par))) par = par[-which(is.na(par))]
   		lko = lk
   		out = lk_comp_latent_cont(Y,yv,Piv,PI,Mu,Si,k)
   		lk = out$lk; Phi = out$Phi; L = out$L; pv = out$pv
# Display output
       	if(it/10 == floor(it/10)){
       		#cat("",sprintf("%11g",c(it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
       		cat(sprintf("%11g",c(k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")
       	}
   		lkv = c(lkv,lk)
	}
	if(it/10 > floor(it/10))  cat(sprintf("%11g",c(k,start,it,lk,lk-lko,max(abs(par-paro)))),"\n",sep=" | ")

# Compute number of parameters
    np = k*r+r*(r+1)/2
    np = np+(k-1)*(nc1+1)
    if(param=="multilogit") np = np+(k-1)*(nc2+1)*k else if(param=="difflogit")  np = np+(k-1)*(nc2+k)
  	aic = -2*lk+np*2
  	bic = -2*lk+np*log(n)
#	out = list(lk=lk,piv=piv,Pi=Pi,Psi=Psi,np=np,aic=aic,bic=bic,lkv=lkv,J=J,V=V1,th=th,sc=sc)
# local decoding
    Ul = matrix(0,n,TT)
    for(i in 1:n) for(t in 1:TT){
    	Ul[i,t] = which.max(V[i,,t])
    }

	Be = matrix(be,nc1+1,k-1)
	if (is.null(nameBe)){
		if(nc1==0) nameBe = c("Intercept") else nameBe = c("intercept",paste("X1",1:nc1,sep=""))
	}else{
		nameBe = c("intercept",nameBe)
	}	

	dimnames(Be) = list(nameBe,logit=2:k)
	if(param=="multilogit"){
		if(is.null(nameGa)){
			if(nc2==0) nameGa = c("Intercept") else nameGa = c("intercept", paste("X2",1:nc2,sep=""))
		}else{
			nameGa = c("intercept",nameGa)
		}
		if(k>2) {
			Ga = array(as.vector(Ga),c(nc2+1,k-1,k))
			dimnames(Ga) = list(nameGa,logit=2:k,logit=1:k)
		}else if(k==2){ 
			
			dimnames(Ga) = 	list(nameGa,logit=1:k)
		}
		
	}else if(param=="difflogit"){
		Ga0 = Ga
		Ga = vector("list",2)
		seGa = vector("list",2)
		Ga[[1]] = t(matrix(Ga0[1:(k*(k-1))],k-1,k))
		Ga[[2]] = matrix(Ga0[(k*(k-1))+(1:((k-1)*nc2))],nc2,k-1)
		if(is.null(nameGa)){
			nameGa2 = paste("X2",1:nc2,sep="")
		}else{
			nameGa2 = nameGa
		}
		if (k==2) {
			dimnames(Ga[[1]]) = list(intercept=1:k,logit=k)
			dimnames(Ga[[2]])=list(nameGa2,logit=k)
		} else if (k>2){
			dimnames(Ga[[1]]) = list(intercept=1:k,logit=2:k)
			dimnames(Ga[[2]])=list(nameGa2,logit=2:k)
		}
	}
	# adjust output
	lk = as.vector(lk)
	if(output){
		dimnames(Piv)=list(subject=1:n,state=1:k)
		dimnames(PI)=list(state=1:k,state=1:k,subject=1:n,time=1:TT)
	}
	if(r==1) dimnames(Mu) = list(item=1,state=1:k) else dimnames(Mu)=list(item=1:r,state=1:k)
	dimnames(Si)=list(item=1:r,item=1:r)
	out = list(lk=lk,Be=Be,Ga=Ga,Mu=Mu,Si=Si,np=np,aic=aic,bic=bic,lkv=lkv,
	           call=match.call(),param=param)
	 
	# final output
	if(miss) out$Y = Y
    if(output){
    	out$PI = PI
    	out$Piv = Piv
    	out$Ul = Ul
    }
    #cat(" |-------------|-------------|-------------|-------------|\n");
    cat("------------|-------------|-------------|-------------|-------------|-------------|\n");
    class(out)="LMlatentcont"
	return(out)
}
