lk_comp_latent_cont <- function(Y,R,Piv,PI,Mu,Si,k,fort=FALSE){

# der = TRUE for derivative
# ---- Preliminaries ----
  if(!fort){
    dnorm1 <-function(y,mu,si2) f = exp(-(y-mu)^2/si2/2)/sqrt(2*pi*si2)
    dmvnorm1 <-function(y,mu,Si) f = exp(-c((y-mu)%*%solve(Si)%*%(y-mu))/2)/sqrt(det(2*pi*Si))
  }
  sY = dim(Y)
  n = as.integer(sY[1])
  TT = as.integer(sY[2])
  if(length(sY)==2) r = 1 else r = sY[3]
  r = as.integer(r)
  if(r==1) if(is.matrix(Y)) Y = array(Y,c(dim(Y),1))
  if(r==1){
    if(is.matrix(Y)) Y = array(Y,c(dim(Y),1))
    if(is.matrix(R)) R = array(R,c(dim(R),1))
  }
  if(is.null(R)) miss = FALSE else miss = any(!R)

  # ---- Compute conditional probabilities ----
  Phi = array(1,c(n,k,TT)); L = array(0,c(n,k,TT))
  if(miss){
    if(fort){
      RR = array(as.integer(1*R),c(n,TT,r))
      out = .Fortran("normmiss",Y,RR,n,TT,r,k,Mu,Si,Phi=Phi)
      Phi = out$Phi
    }else{
      for(i in 1:n) for(t in 1:TT){
        if(all(!R[i,t,])){
          Phi[i,,t] = 1
        }else{
          indo = R[i,t,]
          if(sum(R[i,t,])==1) for(u in 1:k) Phi[i,u,t] = pmax(dnorm1(Y[i,t,][indo],Mu[indo,u],Si[indo,indo]),0.1^300)
          else for(u in 1:k) Phi[i,u,t] = pmax(dmvnorm1(Y[i,t,][indo],Mu[indo,u],Si[indo,indo]),0.1^300)
        }
      }
    }
  }else{
    for(u in 1:k) for(t in 1:TT) Phi[,u,t] =  pmax(dmvnorm(matrix(Y[,t,],n,r),Mu[,u],Si),0.1^300)
  }

# ---- forward recursion ----
#	Phi = pmax(Phi,10^-30)
  L[,,1] = Phi[,,1]*Piv
  for(t in 2:TT){
    # for(u in 1:k) Phi[,u,t] = dmvnorm(matrix(Y[,t,],n,r),Mu[,u],Si)
#		Phi = pmax(Phi,10^-30)
    for(i in 1:n)	L[i,,t] = L[i,,t-1]%*%PI[,,i,t]
    L[,,t] = L[,,t]*Phi[,,t]
  }
  if(n==1) pv = sum(L[,,TT]) else pv = rowSums(L[,,TT])
  lk = sum(log(pv))

# ---- output ----
  out = list(lk=lk,Phi=Phi,L=L,pv=pv)
}
