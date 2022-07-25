#!/usr/bin/env Rscript

path<-getwd()
setwd(path)

# Error tolerance of noise exponent change (used in the original fit).
ef=1.2

# Error tolerance of of A/f^B (from fitting) to the data (used in the new fit).
ef2=0.03

# First and last temperature indices.
i1=1
i2=40

# Sometimes the fit program will have trouble finding the fit range at a given temperature.  The program may exit with an error.  The program may instead get stuck finding the range (it should take 5 seconds at most at a given temperature).  In this case, the program has to be manually terminated.  The index where the program has trouble should be entered into "previous."  It will use the fit range from the previous temperature index.
previous=c(0)
#previous=c(3,7,8,9)

# If the noise exponent is too small (~0.2), the program will have trouble finding the fit range.  Placing those temperature indices into "skip" will skip those temperatures.
skip=c(1:10)
#skip=c(1:6)


# The temperature file is read.  It has a length of 40.
temperature=read.table("Data/1/temperature-1.txt",header=F,col.names=c("T"))

# The first generated frequency fit range of the power spectra as a function of temperature.
sink('original-range-x.dat')
sink()

# The second generated frequency fit range of the power spectra as a function of temperature.
sink('new-range-x.dat')
sink()

# The first fit parameters (amplitude and exponent) as a function of temperature.
sink('original-fit-x.dat')
sink()

# The second fit parameters (amplitude and exponent) as a function of temperature.
sink('new-fit-x.dat')
sink()

# The modified range file used in refit.R
sink('modified-range-x.dat')
sink()

for (k in setdiff(i1:i2,skip))
{
	cat("\n",k)
	
	if (is.element(k,previous)==TRUE)
	{
		r1n=r1p
		r2n=r2p
		cat("\tUsing previous range.")
	}
	else
	{
		Ld3=0
		et=1.0
		r1t=0
		r2t=0
		r1n=0
		r2n=0

		while(Ld3==0)
		{
			name=sprintf("Data/FFT/fftt-10-%d.txt",k)
			datar=read.table(name,header=FALSE, col.names=c("fd","xd","yd","zd"))
			datarx=datar[,1:2]
			L=length(datarx[,1])
			x1=datarx[1,1]
			y1=datarx[1,2]
			x2=datarx[L,1]
			y2=datarx[L,2]
			Ap=100;
			Bp=100;
			Bg=-(log10(y2)-log10(y1))/(log10(x2)-log10(x1))
			Ag=10**(log10(y2)+Bg*log10(x2))

			sink('range-error-x.dat')
			sink()

# Limits on the noise power fit region with initial power (Pi) and final power (Pf):
#	Lower limit  10^Pi
#	Upper limit  10^log10(0.5)=0.5
			Pi=floor(1/log(10,datar$fd[1])+1)
			Pf=log10(0.5)

			Ld=10
			N=2

# A/f^B is fit over equal ranges on a log-log plot.  The percent change (error1) between the previous (Bp) and current exponents (Bn) are calculated.  If the 1/f region, the noise exponent (B) should be relatively constant, so the percent change should be smaller
			Lf=floor((Pf-Pi)*Ld-N+2)
			for (j in seq(from=1,to=((Pf-Pi)*Ld-N+2),by=1)){

				P1=10^(Pi+1.0/Ld*(j-1))
				P2=10^(Pi+1.0/Ld*(j-2+N))

				datarx1=datarx[datarx[,1]<=P2 & datarx[,1]>=P1,]
				colnames(datarx1)<-c("fd","xd")
				Ln=length(datarx1[,1])
				x1n=datarx1[1,1]
				y1n=datarx1[1,2]
				x2n=datarx1[Ln,1]
				y2n=datarx1[Ln,2]
				Bg=-(log10(y2n)-log10(y1n))/(log10(x2n)-log10(x1n))
				Ag=10**(log10(y2n)+Bg*log10(x2n))

				sink('range-error-x.dat',append=TRUE)
				if (Ln>2){
						md<-nls(datarx1$xd~A*datarx1$fd^(-B),data=datarx1,start=list(A=Ag,B=Bg),control=nls.control(maxiter=1000000,minFactor=1e-21))
						An=summary(md)$coefficients[1,1]
						Bn=summary(md)$coefficients[2,1]
				}
				error1=(Bn-Bp)/Bp*100
				if (j==1){
					cat(P1,"\t",P2,"\t",An,"\t",Bn,"\t",sep="")
				} else {
					cat(error1,"\n",P1,"\t",P2,"\t",An,"\t",Bn,"\t",sep="")
				}
				sink()
				Bp=Bn
				Ap=An
			}
			sink('range-error-x.dat',append=TRUE)
				cat(-100)
			sink()

			name2=sprintf("range-error-x.dat")
			data2=read.table(name2,header=FALSE, col.names=c("r1","r2","Ag2","Bg2","error2"))
			data3=data2[abs(data2[,5])<et & data2[,4]>0.2 & data2[,1]>1e-4,]
			Ld3=length(data3[,1])
			r1f=data3$r1[1]
			r2f=data3$r2[length(data3[,1])]
			et=et*ef
		}
		
		if (Ld3==1)
		{
			r1n=r1f
			r2n=r2f
		}
		else
		{
			for (m in seq(from=1,to=Ld3-1,by=1)){
				if (data3[m,2]==data3[m+1,1] & r1n==0){
					r1n=data3[m,1]
				}
				if (data3[m,2]==data3[m+1,1] & r1n>0){
					r2n=data3[m+1,2]
				}
			}
		}

		if (r1n==0 & r2n==0){
			r1n=r1p
			r2n=r2p
		}

	}

	sink('original-range-x.dat',append=TRUE)
		cat(temperature$T[k],"\t",k,"\t",r1n,"\t",r2n,"\n",sep="")
	sink()

	r1p=r1n
	r2p=r2n

	fftr=read.table(sprintf("Data/FFT/fftt-10-%d.txt",k),header=FALSE,col.names=c("fd","xd","yd","zd"))
	fftrn=fftr[fftr[,1]>=r1n & fftr[,1]<=r2n,c(1,2)]
	Ln=length(fftrn[,1])
	x1n=fftrn[1,1]
	y1n=fftrn[1,2]
	x2n=fftrn[Ln,1]
	y2n=fftrn[Ln,2]
	Bg=-(log10(y2n)-log10(y1n))/(log10(x2n)-log10(x1n))
	Ag=10**(log10(y2n)+Bn*log10(x2n))
	md<-nls(fftrn$xd~A*fftrn$fd^(-B),data=fftrn,start=list(A=Ag,B=Bg))
	An=summary(md)$coefficients[1,1]
	Bn=summary(md)$coefficients[2,1]

	sink('original-fit-x.dat',append=T)
		cat(temperature$T[k],"\t",An,"\t",Bn,"\n",sep="")
	sink()


# The frequency range is recalculated based on the error between A/f^B and the data.
	fftrn2=fftr[abs((fftr[,2]/(An/fftr[,1]^Bn)-1)*100)<ef2,c(1,2)]
	r1n=fftrn2[1,1]
	r2n=fftrn2[length(fftrn2[,1]),1]
	if (is.element(k,previous)==TRUE){
		r1n=r1p
		r2n=r2p
	} else if (is.element(k,previous)==FALSE){

# A/f^B will fit well in the linear region on a log-log plot.  Just outside of this region where the power spectra starts to curve, the error between A/f^B and the data will still be low, so the curved portions of the power spectra lie within the fit range.  To account for this, the fit range is reduced
		r1n=r1n*1.3
		r2n=r2n*0.8

# If the lower frequency edge of the fit range ended up being higher in frequency than the higher frequency edge, the lower frequency is reduced, and the higher frequency is increased.
		while(r1n>r2n){
			r1n=r1n*1.2/1.3
			r2n=r2n*0.9/0.8
		}
	}

# If the frequency range is too small, it is set to be the frequency range of the previous temperature.
	if ((abs(r1n-r2n)/(r1n+r2n)*200)<10){
		r1n=r1p
		r2n=r2p
	}

	sink('new-range-x.dat',append=TRUE)
		cat(temperature$T[k],"\t",k,"\t",r1n,"\t",r2n,"\n",sep="")
	sink()

	sink('modified-range-x.dat',append=TRUE)
		cat(temperature$T[k],"\t",k,"\t",r1n,"\t",r2n,"\n",sep="")
	sink()

	fftr=read.table(sprintf("Data/FFT/fftt-10-%d.txt",k),header=FALSE,col.names=c("fd","xd","yd","zd"))
	fftrn=fftr[fftr[,1]>=r1n & fftr[,1]<=r2n,]
	Ln=length(fftrn[,1])
	f1n=fftrn[1,1]
	x1n=fftrn[1,2]
	f2n=fftrn[Ln,1]
	x2n=fftrn[Ln,2]
	Bgx=-(log10(x2n)-log10(x1n))/(log10(f2n)-log10(f1n))
	Agx=10**(log10(x2n)+Bgx*log10(f2n))
	if (is.element(k,skip)==TRUE)
	{
		Agx=Agxp
		Bgx=Bgxp
		Agy=Agyp
		Bgy=Bgyp
		Agz=Agzp
		Bgz=Bgzp
	}
	mdx<-nls(fftrn$xd~Ax*fftrn$fd^(-Bx),data=fftrn,start=list(Ax=Agx,Bx=Bgx),control=nls.control(maxiter=100000,minFactor=1e-21))
	Axn=summary(mdx)$coefficients[1,1]
	Bxn=summary(mdx)$coefficients[2,1]
	sink('new-fit-x.dat',append=T)
		cat(temperature$T[k],"\t",Axn,"\t",Bxn,"\n",sep="")
	sink()
	Agxp=Agx
	Bgxp=Bgx
}

cat("\n")

