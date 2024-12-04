ViewStruc<-function(base.dir,finlist,SubDB,dbkey){
	library(rcdk)
	TMP<-SelectSubDB(base.dir,finlist,SubDB)
	subDB<-TMP[[1]]
	dbkey_row<-which(as.integer(subDB[[1]][,3])%in%as.integer(dbkey))
	if(subDB[[1]][dbkey_row,6]=="NULL"){
	 print("No structure available")
	}else{
	 smile<-as.character(c(subDB[[1]][dbkey_row,6]))
	 mol<-parse.smiles(smile)
	 dep<-get.depictor(width=300,height=300,style="nob")
#	 copy.image.to.clipboard(mol[[1]], dep)
	 img<-view.image.2d(mol[[1]],dep)
	 plot.new()
	 plot.window(xlim=c(0,10),ylim=c(0,10))
	 rasterImage(img,0,0,10,10)
	}
}

# ViewStruc(base.dir,finlist,SubDB="FTneg",dbkey)