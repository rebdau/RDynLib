Local.net<-function(net.lst,dbkey,nettype,nr_of_seq=NULL){
	if (is.null(nr_of_seq)) nr_of_seq=2
	nodes.net<-net.lst[[1]]
	edges.net<-net.lst[[2]]
	sequ=dbkey
	sequ.cont<-as.integer()
	counter=1
	repeat{
	 sequ.row<-which(edges.net[,1]%in%as.integer(sequ))
	 sequ.cont<-c(sequ.cont,sequ.row)
	 sequ<-edges.net[sequ.row,2]
	 if(counter==nr_of_seq)break
	 counter=counter+1
	}
	edges.netloc<-edges.net[sequ.cont,]
	if(dim(edges.netloc)[1]==0){
	 print("no conversions")
	 return()
      }
	nodes.loc<-sort(as.integer(c(edges.netloc[,1],edges.netloc[,2])))
	nodes.loc<-nodes.loc[!duplicated(nodes.loc)]
	nodes.netloc<-nodes.net[nodes.net[,1]%in%nodes.loc,]
	library(RColorBrewer)
	pal3<-brewer.pal(8,"BuGn")
	library(igraph)
	net<-graph.data.frame(edges.netloc,nodes.netloc,directed=T)	
	if (nettype=="CSPP") {
	 plot(net,edge.arrow.size=.6,vertex.label.cex=0.7,vertex.size=V(net)$X2,
		edge.label=E(net)$conv.type,edge.label.cex=0.7,
		edge.label.color=1,edge.color=pal3[E(net)$dot.rel],
		edge.width=40*E(net)$common.nr.rel,main="CSPP")
	} else if (nettype=="GNPS") {
	 plot(net,edge.arrow.size=.6,vertex.label.cex=0.7,vertex.size=V(net)$X2,
		edge.label=E(net)$mass.diff,edge.label.cex=0.7,
		edge.label.color=1,edge.color=pal3[E(net)$dot.rel],
		edge.width=40*E(net)$common.nr.rel,main="GNPS-like")
	} else {
	 print("Define type of net: CSPP or GNPS")
	}
	return(list(nodes.netloc,edges.netloc))
}

# locNET.list<-Local.net(net.lst,dbkey=8,nettype="CSPP")
# locNET.list<-Local.net(net.lst,dbkey=8,nettype="GNPS")
