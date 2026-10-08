mc_multinom <- function(formula,data,baseline=c("initial","central"),
                        model_int=c("all","const","dist1","dist2","dist3","two","symm","rsymm","diff"),
                        model_cov=c("all","const","dist1","dist2","dist3","two","symm","rsymm","diff"),
                        formula_init=NULL, reltol=10^-10,output=FALSE){

#---- preliminaries ----
  model_int = match.arg(model_int)
  model_cov = match.arg(model_cov)
  baseline = match.arg(baseline)
  XX = model.matrix(as.formula(formula),data)
  if(is.null(formula_init)){
    XX1 = XX
  }else{
    XX1 = model.matrix(as.formula(formula_init),data)
  }
  yname = all.vars(formula)[1]
  yv = data[,names(data)==yname]
  k = max(yv)
  idv = data$id
  tv = data$t
  id1v = idv[tv==1]
  XX1 = XX1[tv==1,,drop=FALSE]
  y1v = yv[tv==1]
  n = length(unique(idv))
  Tv = tapply(tv,idv,max)
  TT = max(Tv)
  nc1 = ncol(XX1)-1
  nc = ncol(XX)-1
  knc1 = (k-1)*(1+nc1)
  knc = (k-1)*(1+nc)
  out = design_matrices(k,nc,baseline,model_int,model_cov)
  G = out$G; Z = out$Z; GG = out$GG

#---- starting values ----
  la = rep(0,knc1)
  psi = rep(0,ncol(Z))
  eta = Z%*%psi

#---- initial probabilities log-likelihood ----
  lk1 = 0
  for(i in 1:n){
    Xi = diag(k-1)%x%t(XX1[i,])
    tmp = exp(c(G%*%Xi%*%la))
    tmp = tmp/sum(tmp)
    lk1 = lk1+log(tmp[y1v[i]])
  }
  cat("------------|-------------|-------------|\n")
  cat("    step    |     lk1     |   lk1-lk1o  |\n")
  cat("------------|-------------|-------------|\n")
  cat(sprintf("%11g", c(0, lk1)), "\n", sep = " | ")
#---- iterate until convergence ----
  lk1o = lk1; it = 0
  while((lk1-lk1o)/max(abs(lk1o),10^-300)>reltol || it==0){
    lk1o = lk1; it = it+1

# score vector and information matrix
    sv1 = rep(0,knc1)
    J1 = matrix(0,knc1,knc1)
    for(i in 1:n){
      ri = rep(0,k); ri[y1v[i]] = 1
      Xi = diag(k-1)%x%t(XX1[i,])
      GXi = G%*%Xi
      tmp = exp(c(GXi%*%la))
      tmp = tmp/sum(tmp)
      sv1 = sv1+t(GXi)%*%(ri-tmp)
      Omi = diag(tmp)-tmp%o%tmp
      J1 = J1+t(GXi)%*%Omi%*%GXi
    }

# update parameters
    dla = try(solve(J1,sv1))
    if(inherits(dla,"try-error")) dla = ginv(J1)%*%sv1
    la = la+dla

# initial probabilities log-likelihood
    lk1 = 0
    for(i in 1:n){
      Xi = diag(k-1)%x%t(XX1[i,])
      tmp = exp(c(G%*%Xi%*%la))
      tmp = tmp/sum(tmp)
      lk1 = lk1+log(tmp[y1v[i]])
    }
    cat(sprintf("%11g", c(it, lk1, lk1-lk1o)), "\n", sep = " | ")
  }
  cat("------------|-------------|-------------|\n")
  
#---- transition probabilities log-likelihood ----
  lk2 = 0
  j = 0
  for(i in 1:n){
    j = j+1
    for(t in 2:Tv[i]){
      j = j+1
      Xit = diag(k-1)%x%t(XX[j,])
      u = yv[j-1]
      ind = (u-1)*knc+(1:knc)
      tmp = exp(c(GG[,,u]%*%Xit%*%eta[ind]))
      tmp = tmp/sum(tmp)
      lk2 = lk2+log(tmp[yv[j]])
    }
  }
  cat("------------|-------------|-------------|\n")
  cat("    step    |     lk2     |   lk2-lk2o  |\n")
  cat("------------|-------------|-------------|\n")
  cat(sprintf("%11g", c(0, lk2)), "\n", sep = " | ")

#---- iterate until convergence ----
  lk2o = lk2; it = 0
  while((lk2-lk2o)/abs(lk2o)>reltol || it==0){
    lk2o = lk2; it = it+1

# score vector and information matrix
    sv2 = rep(0,k*knc)
    J2 = matrix(0,k*knc,k*knc)
    j = 0
    for(i in 1:n){
      j = j+1
      for(t in 2:Tv[i]){
        j = j+1
        ri = rep(0,k); ri[yv[j]] = 1
        Xit = diag(k-1)%x%t(XX[j,])
        u = yv[j-1]
        ind = (u-1)*knc+(1:knc)
        GXit = GG[,,u]%*%Xit
        tmp = exp(c(GXit%*%eta[ind]))
        tmp = tmp/sum(tmp)
        sv2[ind] = sv2[ind]+t(GXit)%*%(ri-tmp)
        Omi = diag(tmp)-tmp%o%tmp
        J2[ind,ind] = J2[ind,ind]+t(GXit)%*%Omi%*%GXit
      }
    }
    sv2 = t(Z)%*%sv2
    J2 = t(Z)%*%J2%*%Z

# update parameters
    dpsi = try(solve(J2,sv2))
    if(inherits(dpsi,"try-error")) dpsi = ginv(J2)%*%sv2
    psi = psi+dpsi
    eta = Z%*%psi
    
# transition probabilities log-likelihood
    lk2 = 0
    j = 0
    for(i in 1:n){
      j = j+1
      for(t in 2:Tv[i]){
        j = j+1
        Xit = diag(k-1)%x%t(XX[j,])
        u = yv[j-1]
        ind = (u-1)*knc+(1:knc)
        tmp = exp(c(GG[,,u]%*%Xit%*%eta[ind]))
        tmp = tmp/sum(tmp)
        lk2 = lk2+log(tmp[yv[j]])
      }
    }
    cat(sprintf("%11g", c(it, lk2, lk2-lk2o)), "\n", sep = " | ")
  }
  cat("------------|-------------|-------------|\n")

#---- compute final probabilities ----
  if(output){
    Piv = matrix(0,n,k)
    for(i in 1:n){
      Xi = diag(k-1)%x%t(XX1[i,])
      tmp = exp(c(G%*%Xi%*%la))
      Piv[i,] = tmp/sum(tmp)
    }
    PI = array(0,c(k,k,n,TT))
    j = 0
    for(i in 1:n){
      j = j+1
      for(t in 2:Tv[i]){
        j = j+1
        Xi = diag(k-1)%x%t(XX[j,])
        for(u in 1:k){
          ind = (u-1)*knc+(1:knc)
          tmp = exp(c(GG[,,u]%*%Xi%*%eta[ind]))
          PI[u,,i,t] = tmp/sum(tmp)
        }
      }
    }
  }

#---- output ----
  lk = lk1+lk2
  np = knc1+ncol(Z)
  aic = -2*lk+2*np
  bic = -2*lk+log(n)*np
  sela = sqrt(diag(solve(J1)))
  iJ2 = try(solve(J2))
  if(inherits(iJ2,"try-error")) iJ2 = ginv(J2)
  if(any(diag(J2)<0)) print("negative elements in the diagonal of the inverse information")
  sepsi = sqrt(abs(diag(iJ2)))
  out = list(la=la,psi=psi,eta=eta,lk=lk,lk1=lk1,lk2=lk2,np=np,aic=aic,bic=bic,
             sela=sela,sepsi=sepsi)
  out$call = match.call()
  if(output) out = c(out,list(Piv=Piv,PI=PI))
  class(out)="mc_multinom"
  return(out)

}
