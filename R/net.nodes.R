net.nodes<-function(row.start,row.end,min,subDB){
	nodes.net<-data.frame(matrix(ncol=2,nrow=row.end-row.start+1))
	i=row.start
	repeat{
	 pMS2i<-as.integer(unlist(strsplit(subDB[[3]][[i]],split=",")))
	 pMS2int<-as.numeric(unlist(strsplit(subDB[[4]][[i]],split=",")))
	 ac<-data.frame(cbind(pMS2i,pMS2int))
	 ac<-ac[!ac[,2]<min,]
	 nodes.net[i,1]<-subDB[[1]][i,3]
	 nodes.net[i,2]<-dim(ac)[1]
	 if(i==row.end)break
	 i=i+1
	}
	nodes.net<-nodes.net[row.start:row.end,]
	return(nodes.net)
}

