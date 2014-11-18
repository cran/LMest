summary.LMlatent<-function(object,...){ 
cat("Call:\n")
print(object$call)
cat("\nCoefficients:\n")
	cat("\n Be - Parameters affecting the logit for the initial probabilities:\n")
	print(object$Be)
	if(is.null(object$seBe)==FALSE){
		cat("\n Standard errors for Be:\n")
		print(object$seBe)
	}
	cat("\n Ga - Parameters affecting the logit for the transition probabilities:\n")
	print(object$Ga)	
	if(is.null(object$seGa)==FALSE){
		cat("\n Standard errors for Ga:\n")
		print(object$seGa)		
	}
	cat("\n Psi - Conditional response probabilities:\n")
	print(object$Psi)
	if(is.null(object$sePsi)==FALSE){
		cat("\n Standard errors for the conditional response probability matrix:\n")
		print(object$sePsi)
	}
}