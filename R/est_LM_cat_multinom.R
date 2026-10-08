est_LM_cat_multinom <- function(resp_name, formula, k, data, reltol = 10^-10,
                                rand.start = FALSE,
                                nstart = 1, baseline = c("initial", "central"),
                                model_int = c("all", "const", "dist1", "dist2",
                                            "dist3", "two", "symm", "rsymm", "diff"),
                                model_cov=c("all", "const", "dist1", "dist2", "dist3",
                                            "two", "symm", "rsymm", "diff"),
                                formula_init = NULL,
                                par.init = list(la = NULL, psi = NULL, Pr = NULL),
                                output = FALSE, fort = TRUE){

#---- preliminaries ----
  silent_try = TRUE

#---- multiple-try ----
  if(nstart>1){
    cat(1,"/",nstart,"\n")
    out = est_LM_cat_multinom(resp_name, formula, k, data, reltol,
                              rand.start = FALSE, nstart = 1, baseline,
                              model_int, model_cov, formula_init, par.init,
                              output, fort)
    lkv = out$lk
    for(it in 2:nstart){
      cat(it,"/",nstart,"\n")
      outr = est_LM_cat_multinom(resp_name, formula, k, data, reltol,
                                 rand.start = TRUE, nstart = 1, baseline,
                                 model_int, model_cov, formula_init, par.init,
                                 output, fort)
      lkv = c(lkv,outr$lk)
      print(sort(lkv))
      if(outr$lk>out$lk) out = outr
    }
    out$lkv = lkv
    return(out)
  }

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
  yv = data[,colnames(data)==resp_name]
  idv = data[,1]; iduv = unique(idv)
  tv = data[,2]; tuv = unique(tv)
  id1v = idv[tv==1]
  XX1 = XX1[tv==1,,drop=FALSE]
  nT = length(yv)
  c = max(yv)
  n = length(iduv)
  Tv = tapply(tv,idv,max)
  TT = max(Tv)
  nc1 = ncol(XX1)-1
  nc = ncol(XX)-1
  knc1 = (k-1)*(1+nc1)
  knc = (k-1)*(1+nc)
  if(k>1){
    out = design_matrices(k,nc,baseline,model_int,model_cov)
    G = out$G; Z = out$Z; GG = out$GG
  }

#---- with only one state ----
  if(k==1){
    Pr = matrix(0,1,c)
    for(y in 1:c) Pr[1,y] = sum(yv==y)/nT
    lk = 0
    for(y in 1:c) lk = lk+sum(yv==y)*log(Pr[1,y])
    cat("-------------|-------------|\n")
    cat("      k      |     lk      |\n")
    cat("-------------|-------------|\n")
    cat(sprintf("%11g", c(k, lk)), "\n", sep = " | ")
    cat("------------|-------------|\n")
  }

#---- with more than one state ----
  if(k>1){

#---- starting values ----
    if(rand.start){
      clv = sample(1:k, nT, replace = TRUE)
    }else{
      out = try(kmeans(yv,k,nstart = 500),silent=TRUE)
      if(inherits(out,"try-error")) out = kmeans(yv+seq(0,0.1,length.out = length(yv)),k,nstart = 500)
      clv = out$cluster
    }
    pr = tapply(yv,clv,mean)
    ind = order(pr)
    CL = matrix(0,nT,k); CL[cbind(1:nT,clv)] = 1
    CL = CL[,ind]; clv = CL%*%(1:k)
    if(is.null(par.init$Pr) | rand.start){
      Pr = matrix(0,k,c)
      for(u in 1:k) for(y in 1:c) Pr[u,y] = mean(yv[clv==u]==y)
      Pr = Pr+1/k/2
      Pr = 1/rowSums(Pr)*Pr
    }else{
      Pr = par.init$Pr
    }
    if(is.null(par.init$la) | rand.start){
      la = rep(0,knc1)
    }else{
      la = par.init$la
    }
    Piv = comp_Piv(n,k,XX1,G,la,fort)
    if(is.null(par.init$psi) | rand.start){
      psi = rep(0,ncol(Z))
    }else{
      psi = par.init$psi
    }
    eta = Z%*%psi
    PI = comp_PI(k,n,Tv,XX,GG,eta,fort)

#---- EM algorithm ----
# initial log-likelihood
    Pc = matrix(0,nT,k)
    for(u in 1:k) Pc[,u] = Pr[u,yv]
    out = forward_Multinom(Pc,Piv,PI,n,k,Tv,fort)
    L = out$L; lkv = out$lkv
    lk = sum(lkv)
    cat("------------|-------------|-------------|-------------|\n")
    cat("      k     |     step    |     lk      |    lk-lko   |\n")
    cat("------------|-------------|-------------|-------------|\n")
    cat(sprintf("%11g", c(k, 0, lk)), "\n", sep = " | ")

# iterate until convergence
    lko = lk; it = 0
    while((lk-lko)/abs(lko)>reltol || it==0){
      lko = lk; it = it+1

#---- E-step ----
      out = backward_Multinom(L,Pc,PI,n,k,Tv,fort)
      W = out$W; W2 = out$W2
      W1 = matrix(0,sum(Tv),k)
      j = 0
      for(i in 1:n) for(t in 1:Tv[i]){
        j = j+1
        W1[j,] = W[i,,t]
      }
      

#---- M-step ----
# update conditional response probabilities
      for(u in 1:k){
        wv = W1[,u]
        for(y in 1:c) Pr[u,y] = sum(wv[yv==y])/sum(wv)
      }

# update initial probabilities
      out = comp_sc_Piv(n,k,XX1,G,W,la,fort)
      sv1 = out$sv1; J1 = out$J1
      dla = try(solve(J1,sv1),silent=silent_try)
      if(inherits(dla,"try-error")) dla = ginv(J1)%*%sv1
      la0 = la
      la = la+dla
      Piv = comp_Piv(n,k,XX1,G,la,fort)

# update transition probabilities
      out = comp_sc_PI(k,n,Tv,XX,GG,W2,eta,fort)
      sv2 = out$sv2; J2 = out$J2
      sv2 = t(Z)%*%sv2
      J2 = t(Z)%*%J2%*%Z
      dpsi = try(solve(J2,sv2),silent=silent_try)
      if(inherits(dpsi,"try-error")) dpsi = ginv(J2)%*%sv2
      psi = psi+dpsi
      eta = Z%*%psi
      PI = comp_PI(k,n,Tv,XX,GG,eta,fort)

# recompute log-likelihood
      Pc = matrix(0,nT,k)
      for(u in 1:k) Pc[,u] = Pr[u,yv]
      out = forward_Multinom(Pc,Piv,PI,n,k,Tv,fort)
      L = out$L; lkv = out$lkv
      lk = sum(lkv)
      if(is.na(lk)) browser()
      if(it%%10==0) cat(sprintf("%11g", c(k, it, lk, lk-lko)), "\n", sep = " | ")
    }
    if(it%%10>0) cat(sprintf("%11g", c(k, it, lk, lk-lko)), "\n", sep = " | ")
    cat("------------|-------------|-------------|-------------|\n")
  }

#---- output ----
  dimnames(Pr) = list(u=1:k,category=1:c)
  np = (c-1)*k+knc1
  if(k>1) np = np+ncol(Z)
  aic = -2*lk+2*np
  bic = -2*lk+log(n)*np
  out = list(Pr=Pr,lk=lk,np=np,aic=aic,bic=bic)
  if(k>1){
    CL = apply(W,c(1,3),which.max)
    out = c(out,list(la=la,psi=psi,eta=eta,W=W,CL=CL))
    if(output) out = c(out,list(Piv=Piv,PI=PI))
  }
  out$call = match.call()
  class(out) = "LM_cat_multinom"
  return(out)

}