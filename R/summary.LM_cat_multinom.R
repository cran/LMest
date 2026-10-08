summary.LM_cat_multinom <- function(object,...){

  if(!is.null(object$call)){
    cat("Call:\n")
    print(object$call)
  }
  cat("\nCoefficients:\n")
  
  cat("\n Pr - Emission probabilities:\n")
  print(round(c(object$Pr),4))
  
  if(!is.null(object$la)){
    cat("\n la - Parameters affecting the logit for the initial probabilities:\n")
    print(round(c(object$la),4))
    if(!is.null(object$sela)){
      cat("\n Standard errors for la:\n")
      print(round(object$sela,4))
      pvaluela <- 2*pnorm(abs(c(object$la/object$sela)),lower.tail=FALSE)
      cat("\n p-values for la:\n")
      print(round(pvaluela,4))
    }
  }
  
  if(!is.null(object$psi)){
    cat("\n psi - Parameters affecting the logit for the transition probabilities:\n")
    print(round(c(object$psi),4))
    if(is.null(object$sepsi)==FALSE){
      cat("\n Standard errors for Ga:\n")
      print(round(object$sepsi,4))
      pvaluepsi <- 2*pnorm(abs(c(object$psi/object$sepsi)),lower.tail=FALSE)
      cat("\n p-values for psi:\n")
      print(round(pvaluepsi,4))
    }
  }
  
}
