bootstrap_lm_cov_latent <- function(X1,X2,param="multilogit",Psi,Be,Ga,B=100,fort=TRUE){
    
  mPsi = 0
  mBe = 0
  mGa = 0
  m2Psi = 0
  m2Be = 0
  m2Ga = 0
  dPsi = dim(Psi)
  if(length(dPsi)==1) k = 1
  else k = dPsi[2]
  for (b in 1:B) {
    cat("boostrap sample n. ",b,"\n")
    out = draw_lm_cov_latent(X1,X2,param,Psi,Be,Ga,fort=fort)
    Yb = out$Y
    out = est_lm_cov_latent(Yb,X1,X2,param=param,k=k,fort=fort)
    mPsi = mPsi + out$Psi/B
    mBe = mBe + out$Be/B
    mGa = mGa + out$Ga/B
    m2Psi = m2Psi + out$Psi^2/B
    m2Be = m2Be + out$Be^2/B
    m2Ga = m2Ga + out$Ga^2/B
  }
  sePsi = sqrt(m2Psi - mPsi^2)
  seBe = sqrt(m2Be - mBe^2)
  seGa = sqrt(m2Ga - mGa^2)
  out = list(mPsi = mPsi, mBe = mBe, mGa = mGa, sePsi = sePsi, 
             seBe = seBe, seGa = seGa)
  
}