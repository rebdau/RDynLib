Plot_LCalign<-function(LCal.FT,LCal.Syn,startpoint){
	rg<-Regression_LCalign(LCal.FT,LCal.Syn,startpoint)
	plot(LCal.FT~LCal.Syn,xlab="First Exp",ylab="Second Exp",mgp=c(2,0.5,0),col=5,cex.axis=0.8)
	min.FT=max(LCal.FT)/1000
	x.p=seq(min.FT,max(LCal.FT),min.FT)
	t1.p=as.numeric(rep(0,1000))
	t2.p=as.numeric(rep(0,1000))
	for (i in 1:1000){
	 if(rg[3]!=0){
	  t1.p[i]=ifelse(x.p[i]>rg[3],x.p[i]-rg[3],0)
	 }
	 if(rg[5]!=0){
	  t2.p[i]=ifelse(x.p[i]>rg[5],x.p[i]-rg[5],0)
	 }
	}
	y.p=rg[1]+rg[2]*x.p+rg[4]*t1.p+rg[6]*t2.p
	lines(y.p,x.p) 
	pst<-paste("IWLS: C13=",round(rg[1],digits=2),"+",round(rg[2],digits=2),"FT+",
		round(rg[4],digits=2),"FT[time>",round(rg[3],digits=0),"]+",
		round(rg[6],digits=2),"FT[time>",round(rg[5],digits=0),"]")
	title(main="First Exp vs Second Exp retention times (min)",sub=pst,cex.main=1,cex.sub=0.8)
	return(rg)
}

