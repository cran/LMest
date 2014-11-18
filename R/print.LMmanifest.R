print.LMmanifest <-function(x,...){
cat("Call:\n")
print(x$call)
cat("\nConvergence info:\n")
print(cbind(LogLik=x$lk,np=x$npa,AIC=x$aic, BIC=x$bic))
}