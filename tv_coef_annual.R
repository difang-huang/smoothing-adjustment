
library(tvReg)

library(readxl)
dguep <- read_excel("dguep_annual.xls")

n0=dim(dguep)[1]
y=dguep$dg1
x=dguep$uep1

beta=tvLM(y~x,bw=0.03,tkernel="Gaussian",est="ll")
beta1=beta$coefficients[,2]
#plot(beta1[121:n0-1-120])


ci=confint(beta,tboot = "wild2",level=0.95,runs=500)
blower=ci$Lower[,2]
bupper=ci$Upper[,2]
se=(bupper-blower)/1.96/2
tstat=beta1/se

timerange=seq(1873,2019,length.out = n0-1)

timerangep=timerange[20:(n0-1)]
beta1p=beta1[10:(n0-1-10)]
blowerp=blower[10:(n0-1-10)]
bupperp=bupper[10:(n0-1-10)]
sep=se[10:(n0-1-10)]
tstatp=tstat[10:(n0-1-10)]

plot(timerangep,beta1p,type="l",ylim=c(-0.2,0.9),xlab = "Year",ylab = "Time-varying coefficient",lwd=2)
lines(timerangep,bupperp,col=2,lty=2,lwd=2)
lines(timerangep,blowerp,col=2,lty=2,lwd=2)
abline(h=0)

plot(timerangep,sep,type="l",xlab = "Year",ylab = "standard error",lwd=2,ylim=c(0,0.2))
abline(h=0,lwd=2,lty=2,col=2)


plot(timerangep,tstatp,type="l",xlab = "Year",ylab = "t-statistic",lwd=2,ylim=c(0,28))
abline(h=1.96,lwd=2,lty=2,col=2)

write.csv(cbind(beta1,bupper,blower,tstat),"results_dg1uep1_annual.csv")




