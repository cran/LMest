recursions <-
function(S,yv,Psi,piv,Pi,k,lth,Am,Bm,Cm,b,mod){
	
# forward recursion
ns = dim(S)[1]
TT = dim(S)[2]
M = array(0,c(k,TT,ns))
for(i in 1:ns) for(t in 1:TT) M[,t,i] = Psi[S[i,t]+1,,1]
Q = array(0,c(k,TT,ns))
pv = rep(0,ns)
for(i in 1:ns){
	q = M[,1,i]*piv
	Q[,1,i] = q
	for(t in 2:TT){
		q = M[,t,i]*(t(Pi[,,t])%*%q)
		Q[,t,i] = q	
	}
	pv[i] = sum(q)	
}

#log-likelihood
lk = sum(yv*log(pv))

# backward recurion
Qb = array(0,c(k,TT,ns))
for(i in 1:ns){
	qb = rep(1,k)
	Qb[,TT,i] = qb
	for(t in seq(TT-1,1,-1)){
		qb = (Pi[,,t+1])%*%(M[,t+1,i]*qb)
		Qb[,t,i] = qb	
	}
}

# posterior probabilities
F1 = Q*Qb
for(i in 1:ns) F1[,,i] = F1[,,i]/pv[i]
F2 = array(0,c(k,k,TT,ns))
QbM = Qb*M
for(i in 1:ns) for(t in 2:TT) F2[,,t,i] = Pi[,,t]*(Q[,t-1,i]%o%QbM[,t,i])/pv[i]

#derivative of the forward recursion
ind = 0
#derivative of Psi
Psid = array(0,c(b+1,k,lth))
for(u in 1:k){
	Om = diag(Psi[,u,1])-Psi[,u,1]%o%Psi[,u,1]
	D = Om%*%Am
	for(j in 1:b){
		ind = ind+1
		Psid[,u,ind] = D[,j]
	}		
}
#derivative of piv
pivd = matrix(0,k,lth)
Om = diag(piv)-piv%o%piv
D = Om%*%Bm
for(j in 1:(k-1)){
	ind = ind+1
	pivd[,ind] = D[,j]
}
#derivative of Pi
if(mod==0){
	Pid = array(0,c(k,k,TT,lth))
	for(t in 2:TT) for(u in 1:k){
		Om = diag(Pi[u,,t])-Pi[u,,t]%o%Pi[u,,t]
		D = Om%*%Cm[,,u]
		for(j in 1:(k-1)){
			ind = ind+1
			Pid[u,,t,ind] = D[,j]
		}	
	}
}
if(mod==1){
	Pid = array(0,c(k,k,lth))
	for(u in 1:k){
		Om = diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2]
		D = Om%*%Cm[,,u]
		for(j in 1:(k-1)){
			ind = ind+1
			Pid[u,,ind] = D[,j]
		}	
	}
}
if(mod>1){
	Pid = array(0,c(k,k,2,lth))
	for(u in 1:k){
		Om = diag(Pi[u,,2])-Pi[u,,2]%o%Pi[u,,2]
		D = Om%*%Cm[,,u]
		for(j in 1:(k-1)){
			ind = ind+1
			Pid[u,,1,ind] = D[,j]
		}	
	}	
	for(u in 1:k){
		Om = diag(Pi[u,,mod+1])-Pi[u,,mod+1]%o%Pi[u,,mod+1]
		D = Om%*%Cm[,,u]
		for(j in 1:(k-1)){
			ind = ind+1
			Pid[u,,2,ind] = D[,j]
		}	
	}	
}

#derivative of the forward recursion
Md = array(0,c(k,TT,ns,lth))
for(i in 1:ns) for(t in 1:TT) Md[,t,i,] = Psid[S[i,t]+1,,]
Qd = array(0,c(k,TT,ns,lth))
pvd = matrix(0,ns,lth)
for(j in 1:lth){
	for(i in 1:ns){
		qd = Md[,1,i,j]*piv+M[,1,i]*pivd[,j]
		Qd[,1,i,j] = qd
		for(t in 2:TT){
			q = Q[,t-1,i]
			if(mod==0) qd = Md[,t,i,j]*(t(Pi[,,t])%*%q)+M[,t,i]*(t(Pid[,,t,j])%*%q)+M[,t,i]*(t(Pi[,,t])%*%qd)
			if(mod==1) qd = Md[,t,i,j]*(t(Pi[,,t])%*%q)+M[,t,i]*(t(Pid[,,j])%*%q)+M[,t,i]*(t(Pi[,,t])%*%qd)
			if(mod>1){
				if(t<=mod) qd = Md[,t,i,j]*(t(Pi[,,t])%*%q)+M[,t,i]*(t(Pid[,,1,j])%*%q)+M[,t,i]*(t(Pi[,,t])%*%qd)
				if(t>mod) qd = Md[,t,i,j]*(t(Pi[,,t])%*%q)+M[,t,i]*(t(Pid[,,2,j])%*%q)+M[,t,i]*(t(Pi[,,t])%*%qd)
			}
			Qd[,t,i,j] = qd	
		}
		pvd[i,j] = sum(qd)	
	}	
}

#score
sc = (yv/pv)%*%pvd

# derivative of backward recurion
Qbd = array(0,c(k,TT,ns,lth))
for(j in 1:lth){
	for(i in 1:ns){
		qbd = rep(0,k)
		Qbd[,TT,i,j] = qbd
		for(t in seq(TT-1,1,-1)){
			qb = Qb[,t+1,i]
			if(mod==0) qbd = (Pid[,,t+1,j])%*%(M[,t+1,i]*qb)+(Pi[,,t+1])%*%(Md[,t+1,i,j]*qb)+(Pi[,,t+1])%*%(M[,t+1,i]*qbd)
			if(mod==1) qbd = (Pid[,,j])%*%(M[,t+1,i]*qb)+(Pi[,,t+1])%*%(Md[,t+1,i,j]*qb)+(Pi[,,t+1])%*%(M[,t+1,i]*qbd)
			if(mod>1){
				if(t>=mod) qbd = (Pid[,,2,j])%*%(M[,t+1,i]*qb)+(Pi[,,t+1])%*%(Md[,t+1,i,j]*qb)+(Pi[,,t+1])%*%(M[,t+1,i]*qbd)
				if(t<mod) qbd = (Pid[,,1,j])%*%(M[,t+1,i]*qb)+(Pi[,,t+1])%*%(Md[,t+1,i,j]*qb)+(Pi[,,t+1])%*%(M[,t+1,i]*qbd)
			}
			Qbd[,t,i,j] = qbd	
		}
	}	
}

# posterior probabilities
F1d = array(0,c(k,TT,ns,lth))
F1t = Q*Qb
for(j in 1:lth){
	F1dj = Qd[,,,j]*Qb+Q*Qbd[,,,j]
	for(i in 1:ns) F1d[,,i,j] = F1dj[,,i]/pv[i]-F1t[,,i]*pvd[i,j]/pv[i]^2
	}
F2d = array(0,c(k,k,TT,ns,lth))
for(j in 1:lth){
	QbMdj = Qbd[,,,j]*M+Qb*Md[,,,j]
	for(i in 1:ns) for(t in 2:TT){
		if(mod==0) F2d[,,t,i,j] = Pid[,,t,j]*(Q[,t-1,i]%o%QbM[,t,i])/pv[i]+Pi[,,t]*(Qd[,t-1,i,j]%o%QbM[,t,i])/pv[i]+
                 	   Pi[,,t]*(Q[,t-1,i]%o%QbMdj[,t,i])/pv[i]-Pi[,,t]*(Q[,t-1,i]%o%QbM[,t,i])*pvd[i,j]/pv[i]^2  
		if(mod==1) F2d[,,t,i,j] = Pid[,,j]*(Q[,t-1,i]%o%QbM[,t,i])/pv[i]+Pi[,,t]*(Qd[,t-1,i,j]%o%QbM[,t,i])/pv[i]+
                 	   Pi[,,t]*(Q[,t-1,i]%o%QbMdj[,t,i])/pv[i]-Pi[,,t]*(Q[,t-1,i]%o%QbM[,t,i])*pvd[i,j]/pv[i]^2                
		if(mod>1){
			if(t <=mod) F2d[,,t,i,j] = Pid[,,1,j]*(Q[,t-1,i]%o%QbM[,t,i])/pv[i]+Pi[,,t]*(Qd[,t-1,i,j]%o%QbM[,t,i])/pv[i]+
                 	   Pi[,,t]*(Q[,t-1,i]%o%QbMdj[,t,i])/pv[i]-Pi[,,t]*(Q[,t-1,i]%o%QbM[,t,i])*pvd[i,j]/pv[i]^2 
            if(t>mod) F2d[,,t,i,j] = Pid[,,2,j]*(Q[,t-1,i]%o%QbM[,t,i])/pv[i]+Pi[,,t]*(Qd[,t-1,i,j]%o%QbM[,t,i])/pv[i]+
                 	   Pi[,,t]*(Q[,t-1,i]%o%QbMdj[,t,i])/pv[i]-Pi[,,t]*(Q[,t-1,i]%o%QbM[,t,i])*pvd[i,j]/pv[i]^2 
		}
	} 
}
out = list(lk=lk,sc=sc,F1=F1,F2=F2,F1d=F1d,F2d=F2d)

}
