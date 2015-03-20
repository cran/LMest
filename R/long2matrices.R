long2matrices <- function(id,X,Y){

# preliminaries
	idu = unique(id)
	n = length(idu)
	Tv = table(id)
	TT = max(Tv)
	X = as.matrix(X)
	nx = ncol(X)
	Y = as.matrix(Y)
	ny = ncol(Y)
# create matrices
	XX = array(NA,c(n,TT,nx))
	YY = array(NA,c(n,TT,ny))
	for(i in 1:n){
		ind = which(id==idu[i])
		for(t in 1:TT){
			XX[i,t,] = X[ind[t],] 
			YY[i,t,] = Y[ind[t],]
		}
	}
# output
	out = list(XX=XX,YY=YY)
	out
  
}