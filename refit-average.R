#!/usr/bin/env Rscript


# This script fits A/f^B to the power spectra using initial guesses from "new-fit.dat" within the fit range "modified-range.dat."  The fit parameters as a function of temperature are written to "refit.dat."


path<-getwd()
setwd(path)

temperature=read.table("Data/1/temperature-1.txt",header=F,col.names=c("T"))

sink('refit.dat')
sink()

i1=1
previousfit=read.table("new-fit.dat",header=FALSE)
previousrange=read.table("modified-range.dat",header=FALSE)

# This sets the number of temperature steps the refit is done on.  If some temperature steps were skipped, this takes that into account.
pr=previousrange[1,2]
i2=length(temperature$T)+1-pr
for (k in seq(from=i1,to=i2,by=1))
{
	cat("\n",k)
	r1n=previousrange[k,3]
	r2n=previousrange[k,4]
	Ag=previousfit[k,2]
	Bg=previousfit[k,3]

	fftr=read.table(sprintf("Data/FFT/fftt-10-%d.txt",k+pr-1),header=FALSE,col.names=c("fd","xd","yd","zd"))
	fftrn=fftr[fftr[,1]>=r1n & fftr[,1]<=r2n,]
	Ln=length(fftrn[,1])
	f1n=fftrn[1,1]
	x1n=fftrn[1,2]
	y1n=fftrn[1,3]
	z1n=fftrn[1,4]
	f2n=fftrn[Ln,1]
	x2n=fftrn[Ln,2]
	y2n=fftrn[Ln,3]
	z2n=fftrn[Ln,4]
	Agx=Ag
	Agy=Ag
	Agz=Ag
	Bgx=Bg
	Bgy=Bg
	Bgz=Bg
	mdx<-nls(fftrn$xd~Ax*fftrn$fd^(-Bx),data=fftrn,start=list(Ax=Agx,Bx=Bgx),control=nls.control(maxiter=100000))
	Axn=summary(mdx)$coefficients[1,1]
	Bxn=summary(mdx)$coefficients[2,1]
	mdy<-nls(fftrn$yd~Ay*fftrn$fd^(-By),data=fftrn,start=list(Ay=Agy,By=Bgy),control=nls.control(maxiter=100000))
	Ayn=summary(mdy)$coefficients[1,1]
	Byn=summary(mdy)$coefficients[2,1]
	mdz<-nls(fftrn$zd~Az*fftrn$fd^(-Bz),data=fftrn,start=list(Az=Agz,Bz=Bgz),control=nls.control(maxiter=100000))
	Azn=summary(mdz)$coefficients[1,1]
	Bzn=summary(mdz)$coefficients[2,1]

	sink('refit.dat',append=T)
		cat(temperature$T[k+pr-1],"\t",(Axn+Ayn+Azn)/3.0,"\t",(Bxn+Byn+Bzn)/3.0,"\n",sep="")
	sink()
	Agxp=Agx
	Bgxp=Bgx
	Agyp=Agy
	Bgyp=Bgy
	Agzp=Agz
	Bgzp=Bgz

}
cat("\n")

