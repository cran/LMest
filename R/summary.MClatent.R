summary.MClatent<-function(object,...){ 
cat("Call:\n")
print(object$call)
cat("\nCoefficients:\n")

cat("\n Be - Parameters affecting the logit for the initial probabilities:\n")
print(round(object$Be,4))
if(is.null(object$seBe)==FALSE){
	cat("\n Standard errors for Be:\n")
	print(round(object$seBe,4))
}

cat("\n Ga - Parameters affecting the logit for the transition probabilities:\n")
print(round(object$Ga,4))	
if(is.null(object$seGa)==FALSE){
	cat("\n Standard errors for Ga:\n")
	print(round(object$seGa,4))		
}	

}