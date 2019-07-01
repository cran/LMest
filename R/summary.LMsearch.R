summary.LMsearch<-function(object,...){
  if(!is.null(object$call))
  {
    cat("Call:\n")
    print(object$call)
  }
	kdim = length(object$kv)
	np = lk = AIC = BIC= rep(0,kdim)
	cont=0
	kv <- object$k
	for(k in kv){
		cont=cont+1
		np[cont] = object$out.single[[k]]$np
		lk[cont] = object$lkv[k]
		#AIC[cont] = object$Aic[k]
		#BIC[cont] = object$Bic[k]
	}
	if(!is.null(object$aicv))
	{
	  Aic = object$aicv
	  Bic = object$bicv
	}else{
	  Aic = object$Aic
	  Bic = object$Bic
	}
	print(cbind(states=object$k,lk=lk,np=np,AIC=Aic,BIC=Bic))

}
