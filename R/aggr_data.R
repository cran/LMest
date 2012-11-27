aggr_data <-
function(data){

    # find distinct response configurations and corresponding frequencies
    data_dis = data
    rim = nrow(data)    
    nc = ncol(data)
    freq = rep(0,rim)
    label = rep(0,rim)
    j = 0
    label0 = 1:rim
   	cat("------------|-------------|\n");
   	cat("  remaining |   average   |\n")
   	cat("------------|-------------|\n");
    while(rim>0){
    	j = j+1
    	data_dis[j,] = data[1,]
		D = t(matrix(data_dis[j,],nc,rim))
	    ind = which(rowSums(D==data)==nc)
	    label[label0[ind]] = j
		freq[j] = length(ind)    
		data = data[-ind,]
		label0 = label0[-ind]
	    rim = length(label0)
	    if(rim==1) data = t(data)
	    cat(sprintf("%11g",c(rim,mean(freq[1:j]))),"\n",sep=" | ")    
    } 
   	cat("------------|-------------|\n");
    data_dis = data_dis[1:j,]
    freq = freq[1:j]
    # output
    out = list(data_dis=data_dis,freq=freq,label=label)
}
