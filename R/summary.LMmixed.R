summary.LMmixed<-function(object,...){ 
cat("Call:\n")
print(object$call)
cat("\nCoefficients:\n")
cat("\nMass probabilities:\n")
print(object$la)
cat("\nInitial probabilities:\n")
print(object$Piv)
cat("\nTransition probabilities:\n")
print(object$Pi)
cat("\nConditional response probabilities:\n")
print(object$Psi)
}