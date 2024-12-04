# m/z feature grouping based on equal retention times and a high Pearson
# correlation across biological replicates.
# y10, y20 and fc are the maximum retention time difference, the minimum
# correlation and a vector of sample groups

FeatureGroupXCMS<-function(dt,y10,y20,fc){
	x=dt					# xcms subset table ordered by tR
	xcms2i=x
	o<-order(xcms2i$name)		# xcms subset table ordered by name
	xcms2<-xcms2i[o,]
	fc1=x[,-(1:13)]			# xcms subset table of only abundances
	fc2=1					# lines
	fc3=max(fc)				# total number of lines
	fc4=1					# column with first rep of line
	repeat{
		fc5<-length(fc[fc==fc2])	
			# fc5: number of reps for a particular line
		fc6<-fc4+fc5-1		# column with last rep of line
		x1=t(fc1[,fc4:fc6])	
			# x1: subset of fc1 that contains a particular line
		mu=colMeans(x1,na.rm=T)
		k=ncol(x1)
		nm=as.character(x$name)
		colnames(x1)=nm
		z2=cor(x1,use="p")
		z2=z2[lower.tri(z2)]
		z1=outer(x$rtmed,x$rtmed,"-")
		z1=abs(z1[lower.tri(z1)])
		jj=(1:(k*(k-1)/2))[(z1<y10)&(!is.na(z1))&(z2>y20)&(!is.na(z2))]
		rw=as.integer(c())
		i=1
		while (i<k){
			rw.i<-seq(i+1,k,1)
			rw<-append(rw,rw.i)
			i=i+1
		}
		cw=as.integer(c())
		i=1
		ki<-k-1
		while (i<k){
			cw.i<-rep(i,ki)
			cw<-append(cw,cw.i)
			i=i+1	
			ki=ki-1
		}
		res1<-array(dim=c(0,2))
		i=1
		while (i<length(jj)+1){
			cord.i<-c(rw[jj[i]],cw[jj[i]])
			res1<-rbind(res1,cord.i)
			i=i+1
		}
		res=res1
		grp=2
		Grp<-as.integer(rep(0,length(nm)))
		nm.Grp<-data.frame(cbind(nm,Grp))
		nm.Grp$Grp<-as.integer(nm.Grp$Grp)
		repeat {
			set=as.integer(res[1,])
			l1=length(set)
			repeat {
		  		set=unique(c(set,res[res[,2]%in%set,1])) 
		  		set=unique(c(set,res[res[,1]%in%set,2]))
		  		l2=length(set)
		  		if (l2 == l1) break
		  		l1=l2
			}
			nm.Grp[set,2]=grp
			res<-res[!(res[,1]%in%set),]
			clsres<-is.matrix(res)
			if (clsres=="FALSE") {
				break
			}else{
				if (clsres=="TRUE"){
					if (nrow(res)<1) break
				}
			}
			grp=grp+1	
		}
		o<-order(nm.Grp$nm)
		nm.Grpo<-nm.Grp[o,]
		xcms2<-cbind(xcms2,nm.Grpo$Grp)
		if (fc2==fc3) break
		fc2=fc2+1
		fc4=fc4+fc5
	}
	nrcr=ncol(xcms2)-fc3
	ii<-do.call(order,data.frame(xcms2[,-(1:nrcr)]))
	xcms2<-xcms2[ii,]
	xcms2.cr<-xcms2[,-(1:nrcr)]
	i=2
	j=1
	l=ncol(xcms2.cr)
	fin=c(rep(NA,l))
	fin.c=c(rep(T,l))
	con=c(rep(1,length(xcms2.cr[,1])))
	repeat{
		k=1
		repeat{
			if (xcms2.cr[i,k]==xcms2.cr[j,k]|xcms2.cr[j,k]==1) {
				fin[k]=T
			}else{
				fin[k]=F
			}
			if (k==l) break
			k=k+1
		}
		if (all(fin==fin.c)) {
			con[i]=con[j]
		}else{
			con[i]=con[j]+1
		}
		if (i==length(xcms2.cr[,1])) break
		i=i+1
		j=j+1
	}
	xcms2<-cbind(xcms2,con)
	return(xcms2)
}
