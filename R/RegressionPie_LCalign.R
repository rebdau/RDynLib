#Regression between FT and Synapt chromatograms obtained for the same
#sample. Deviations from a straight curve might occur as different LC
#systems are involved. Therefore, a piecewise regression is performed
#including at most two knots. A hyperbolic curve was noted when
#aligning a Thermo Accela chromatogram with a Waters UPLC chromatogram.
#By looking at the result from multiple experiments, robust regression
#was more accurate than traditional least squares. The 'startpoint'
#allows to decide on the FT retention time from which regression
#should be started. 
 
RegressionPie_LCalign<-function(LCal,startpoint){
	knottime1<-c()
	knottime2<-c()
	VarExp<-c()
	y=as.numeric(LCal[,4]) # Synapt retention times
	x=as.numeric(LCal[,2]) # FT retention times
	out<-lm(y~x)
	knottime1<-append(knottime1,0)
	knottime2<-append(knottime2,0)
	VarExp<-append(VarExp,summary(out)$adj.r.squared)
	i=startpoint
	while(i<as.integer(max(x))){
	 t1<-ifelse(x>i,x-i,0)
	 out<-lm(y~x+t1)
	 knottime1<-append(knottime1,i)
	 knottime2<-append(knottime2,0)
	 VarExp<-append(VarExp,summary(out)$adj.r.squared)
	 i=i+1
	}
	i=2
	while(i<as.integer(max(x))){
	 t1<-ifelse(x>i,x-i,0)
	 j=i+1
	 while(j<as.integer(max(x))){
	  t2<-ifelse(x>j,x-j,0)
	  out<-lm(y~x+t1+t2)
	  knottime1<-append(knottime1,i)
	  knottime2<-append(knottime2,j)
	  VarExp<-append(VarExp,summary(out)$adj.r.squared)
	  j=j+1
	 }
	 i=i+1
	}
	best.model<-data.frame(knottime1,knottime2,VarExp)
	best.model<-best.model[order(best.model[,3],decreasing=T),]
	t1<-ifelse(x>best.model[1,1],x-best.model[1,1],0)
	t2<-ifelse(x>best.model[1,2],x-best.model[1,2],0)
	library(MASS)
	if((best.model[1,1]==0)&(best.model[1,2]==0)){
	 exp.t1=0
	 exp.t2=0
	 out<-rlm(y~x)
	 exp.intercept<-summary(out)$coef[1,1]
	 exp.slope<-summary(out)$coef[2,1]
	}else{
	 if(best.model[1,2]==0){
	  exp.t2=0
	  out<-rlm(y~x+t1)
	  exp.intercept<-summary(out)$coef[1,1]
	  exp.slope<-summary(out)$coef[2,1]
	  exp.t1<-summary(out)$coef[3,1]
	 }else{
	  out<-rlm(y~x+t1+t2)
	  exp.intercept<-summary(out)$coef[1,1]
	  exp.slope<-summary(out)$coef[2,1]
	  exp.t1<-summary(out)$coef[3,1]
	  exp.t2<-summary(out)$coef[4,1]
	 }
	}
	rg<-c(exp.intercept,exp.slope,best.model[1,1],exp.t1,best.model[1,2],exp.t2)
	return(rg)
}

