summary.LMbasic<-function(object,...){ 
piv = cbind(Est_piv = object$piv)
cat("Call:\n")
print(object$call)
cat("\nCoefficients:\n")
cat("\nInitial probabilities:\n")
print(round(piv,4))
if(is.null(object$sepiv)==FALSE){
	cat("\nStandard errors for the initial probabilities:\n")	
	print(round(cbind(se_piv = object$sepiv),4))
}
TT = dim(object$Pi)[3]
cat("\nTransition probabilities:\n")
print(round(object$Pi[,,2:TT],4))
if(is.null(object$sePi)==FALSE){
	cat("\nStandard errors for the transition probabilities:\n")
	print(round(object$sePi[,,2:TT],4))
}	
cat("\nConditional Response probabilities:\n")
print(round(object$Psi,4))
if(is.null(object$sePsi)==FALSE){
	cat("\nStandard errors for the conditional Response probabilities:\n")
	print(round(object$sePsi,4))
}
}