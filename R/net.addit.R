net.addit<-function(base.dir,finlist,SubDB,exp.id,nettype,nr_col=NULL,
				min=NULL,thr1=NULL,thr2=NULL,thr3=NULL){
	if (is.null(thr1)) thr1=3
	if (is.null(thr2)) thr2=0.1
	if (is.null(thr3)) thr3=0.2
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
	compid.start=subDB[[1]][subDB[[1]][,4]==exp.id,3][1]
	row.start=which(subDB[[1]][,3]==compid.start)
	compid.end=subDB[[1]][subDB[[1]][,4]==exp.id,3][length(subDB[[1]][subDB[[1]][,4]==exp.id,3])]
	row.end=which(subDB[[1]][,3]==compid.end)
	if (is.null(min)) {
	 if (TMP[[2]]=="FT"){
	  min=100
	 } else if (TMP[[2]]=="QTOF"){
	  min=5
	 } else {
	  print("Unknown MS Analyzer, set SubDB equal to FTneg,
					 FTpos or QTOFneg")
	 }
	}
	nodes.net<-net.nodes(row.start,row.end,min,subDB)
	setwd(TMP[[4]])
	if (nettype=="CSPP") {
	 if (is.null(nr_col)) nr_col=35
	 comp_add<-read.table("compound_add.txt",header=T,sep="\t",
					stringsAsFactors=FALSE)
	 edges.net<-allCSPP(row.start=compid.start,row.end=compid.end,nr_col,comp_add,thr1=thr1,thr2=thr2,thr3=thr3)
	} else if (nettype=="GNPS") {
	 if (is.null(nr_col)) nr_col=6
	 gnps_add<-read.table("gnps_add.txt",header=T,sep="\t",
					stringsAsFactors=FALSE)
	 edges.net<-allGNPS(row.start=compid.start,row.end=compid.end,nr_col,gnps_add,thr1=thr1,thr2=thr2,thr3=thr3)
	} else {
	 print("Type of net should be CSPP or GNPS.")
	}	
	common.nr<-edges.net[,4]*edges.net[,5]
	edges.net<-data.frame(edges.net,common.nr)
	common.nr.rel<-edges.net[,7]/max(edges.net[,7])
	edges.net<-data.frame(edges.net,common.nr.rel)
	dot.rel<-round(edges.net[,6]*8,digits=0)
	edges.net<-data.frame(edges.net,dot.rel)
	return(list(nodes.net,edges.net))
}

# net.lst<-net.addit(base.dir,finlist,SubDB="FTneg",exp.id=1,nettype="CSPP")
# net.lst<-net.addit(base.dir,finlist,SubDB="FTneg",exp.id=1,nettype="GNPS")


