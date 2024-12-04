Overall.net<-function(base.dir,finlist,SubDB,exp.id,dbkey,nr_of_seq=NULL,thr1=NULL,thr2=NULL,thr3=NULL){
	if (is.null(thr1)) thr1=3
	if (is.null(thr2)) thr2=0.1
	if (is.null(thr3)) thr3=0.2
	if (is.null(nr_of_seq)) nr_of_seq=2	
	oldpar<-par(no.readonly=T)
	par(mfrow=c(1,2),mar=c(1,1,1,1))
	net.lst<-net.addit(base.dir,finlist,SubDB,exp.id,nettype="CSPP",thr1=thr1,thr2=thr2,thr3=thr3)
	locNET1.list<-Local.net(net.lst,dbkey,nettype="CSPP",
						nr_of_seq=nr_of_seq)
	net.lst<-net.addit(base.dir,finlist,SubDB,exp.id,nettype="GNPS",thr1=thr1,thr2=thr2,thr3=thr3)
	locNET2.list<-Local.net(net.lst,dbkey,nettype="GNPS",
						nr_of_seq=nr_of_seq)
	par(oldpar)
	return(list(locNET1.list,locNET2.list))
}

# locNET.list<-Overall.net(base.dir,finlist,SubDB="FTneg",exp.id=1,dbkey=8,nr_of_seq=1)