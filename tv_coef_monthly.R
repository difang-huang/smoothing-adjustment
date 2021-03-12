
library(tvReg)

library(readxl)
dguep <- read_excel("dguep_month.xls")

n0=dim(dguep)[1]
y=dguep$dg2
x=dguep$uep2

beta=tvLM(y~x,bw=0.03,tkernel="Gaussian",est="ll")
beta1=beta$coefficients[,2]
#plot(beta1[121:n0-1-120])


ci=confint(beta,tboot = "wild2",level=0.95,runs=250)
blower=ci$Lower[,2]
bupper=ci$Upper[,2]
se=(bupper-blower)/1.96/2
tstat=beta1/se

timerange=seq(1872,2019,length.out = n0-1)

timerangep=timerange[241:(n0-1)]
beta1p=beta1[121:(n0-1-120)]
blowerp=blower[121:(n0-1-120)]
bupperp=bupper[121:(n0-1-120)]
sep=se[121:(n0-1-120)]
tstatp=tstat[121:(n0-1-120)]

plot(timerangep,beta1p,type="l",ylim=c(0.0023,0.075),xlab = "Year",ylab = "Time-varying coefficient",lwd=2)
lines(timerangep,bupperp,col=2,lty=2,lwd=2)
lines(timerangep,blowerp,col=2,lty=2,lwd=2)
abline(h=0)

plot(timerangep,sep,type="l",xlab = "Year",ylab = "standard error",lwd=2,ylim=c(0,0.01))
abline(h=0,lwd=2,lty=2,col=2)


plot(timerangep,tstatp,type="l",xlab = "Year",ylab = "t-statistic",lwd=2,ylim=c(0,28))
abline(h=1.96,lwd=2,lty=2,col=2)

write.csv(cbind(beta1,bupper,blower,tstat),"results_dg2uep2.csv")




