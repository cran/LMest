design_matrices <- function(k,nc,baseline,model_int,model_cov){

  knc = (k-1)*(1+nc)
  G = diag(k)[,-1]
  GG = array(0,c(k,k-1,k))
  for(u in 1:k){
    if(baseline=="initial") GG[,,u] = G
    if(baseline=="central") GG[,,u] = diag(k)[,-u]
  }
  if(model_int=="all") Z1 = diag(k*(k-1))
  if(model_int=="const") Z1 = matrix(1,k*(k-1),1)
  if(model_int=="dist1"){
    Z1 = matrix(0,k*(k-1),2*(k-1))
    j = 0
    for(u in 1:k) for(v in (1:k)[-u]){
      j = j+1
      Z1[j,k+v-u-1*(v>u)]=1
    }
  }
  if(model_int=="dist2" || model_int=="dist3"){
    Z1 = matrix(0,k*(k-1),k-1)
    j = 0
    for(u in 1:k) for(v in (1:k)[-u]){
      j = j+1
      if(model_int=="dist2") Z1[j,abs(v-u)]=1
      if(model_int=="dist3") Z1[j,abs(v-u)]=2*(v>u)-1
    }
  }
  if(model_int=="two"){
    Z1 = matrix(0,k*(k-1),2)
    j = 0
    for(u in 1:k) for(v in (1:k)[-u]){
      j = j+1
      Z1[j,1+1*(v>u)]=1
    }
  }
  if(model_int=="symm" | model_int=="rsymm"){
    IND = matrix(0,k*(k-1)/2,2)
    j = 0
    for(u in 1:(k-1)) for(v in (u+1):k){
      j = j+1
      IND[j,] = c(u,v)
    }
    Z1 = matrix(0,k*(k-1),k*(k-1)/2)
    j = 0
    for(u in 1:k) for(v in (1:k)[-u]){
      j = j+1
      if(model_int=="symm") Z1[j,which(IND[,1]==min(u,v) & IND[,2]==max(u,v))]=1
      if(model_int=="rsymm") Z1[j,which(IND[,1]==min(u,v) & IND[,2]==max(u,v))]=2*(v>u)-1
    }
  }
  if(model_int=="diff"){
    Z1 = matrix(0,k*(k-1),k-1)
    j = 0
    for(u in 1:k) for(v in (1:k)[-u]){
      j = j+1
      if(v>1) Z1[j,v-1]=1
      if(u>1) Z1[j,u-1]=-1
    }
  }
  if(nc==0) Z = Z1
  if(nc>0){
    if(model_cov=="all") Z2 = diag(k*(k-1))
    if(model_cov=="const") Z2 = rep(1,k*(k-1))
    if(model_cov=="dist1"){
      Z2 = matrix(0,k*(k-1),2*(k-1))
      j = 0
      for(u in 1:k) for(v in (1:k)[-u]){
        j = j+1
        Z2[j,k+v-u-1*(v>u)]=1
      }
    }
    if(model_cov=="dist2" || model_cov=="dist3"){
      Z2 = matrix(0,k*(k-1),k-1)
      j = 0
      for(u in 1:k) for(v in (1:k)[-u]){
        j = j+1
        if(model_cov=="dist2") Z2[j,abs(v-u)]=1
        if(model_cov=="dist3") Z2[j,abs(v-u)]=2*(v>u)-1
      }
    }
    if(model_cov=="two"){
      Z2 = matrix(0,k*(k-1),2)
      j = 0
      for(u in 1:k) for(v in (1:k)[-u]){
        j = j+1
        Z2[j,1+1*(v>u)]=1
      }
    }
    if(model_cov=="symm" || model_cov=="rsymm"){
      IND = matrix(0,k*(k-1)/2,2)
      j = 0
      for(u in 1:(k-1)) for(v in (u+1):k){
        j = j+1
        IND[j,] = c(u,v)
      }
      Z2 = matrix(0,k*(k-1),k*(k-1)/2)
      j = 0
      for(u in 1:k) for(v in (1:k)[-u]){
        j = j+1
        if(model_cov=="symm") Z2[j,which(IND[,1]==min(u,v) & IND[,2]==max(u,v))]=1
        if(model_cov=="rsymm") Z2[j,which(IND[,1]==min(u,v) & IND[,2]==max(u,v))]=2*(v>u)-1
      }
    }
    if(model_cov=="diff"){
      Z2 = matrix(0,k*(k-1),k-1)
      j = 0
      for(u in 1:k) for(v in (1:k)[-u]){
        j = j+1
        if(v>1) Z2[j,v-1]=1
        if(u>1) Z2[j,u-1]=-1
      }
    }
    Z = cbind(Z1%x%diag(1+nc)[,1],Z2%x%diag(1+nc)[,-1])
  }

#---- output ----
  out = list(G=G,Z=Z,GG=GG,Z1=Z1)
  if(nc>0) out$Z2=Z2
  return(out)

}