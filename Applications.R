### ======================================================================= ###
###                             APPLICATION                                 ###  
### ----------------------------------------------------------------------- ###
### Estimation and diagnostic for skew-normal partially linear models       ###
###                                                                         ###
###  Authors: Clécio S. Ferreira and Gilberto A. Paula                      ###
### ======================================================================= ###



################################################################################
# 1. Artificial data
rm(list=ls(all=TRUE))
source('functions.r')

n=200
x=matrix(runif(n),n,1)# parametric explanatory variable
t=seq(0,3*pi,length=n)
ft=matrix(cos(t),n,1)
t=seq(0.6,1.6,length=n)
ft=matrix(cos(4*pi*t)*exp(-t^2/2),n,1)
t0=t
beta0=5
X=x
p=1
sigma2=1
lambda=3
mu=X*beta0+ft

Nc=N(t,t0)
Rc=R(t0)
Qc=Q(t0)
K=Qc%*%solve(Rc)%*%t(Qc)
erro=geraSN(n,0,sigma2,lambda)
y=mu+erro
teta=EMSNspSIC_DF(y,X,Nc,K) #[beta;sigma2;lambda;f;alpha]
mi=miSNsp(teta,y,X,Nc,K) # (beta,sigma2,lambda,f)
mi1=diag(solve(mi))
ep=sqrt(diag(solve(mi1)))

yest=x*as.numeric(teta[1])
vies=sqrt(2/pi*teta[2])*teta[3]/sqrt(1+teta[3]^2)

plot(x,y,type='p')
lines(x,beta0*x,lwd=3,xlab='x',ylab='y')
lines(x,yest,lty=2,lwd=3,col='gray')

fnp=Nc%*%teta[4:(length(teta)-1)] +vies
errof=mi1[4:length(mi1)]
res=y-yest

plot(t,res,type='p')
lines(t,ft+vies,lwd=3,xlab='t',ylab='f(t)')
lines(t,fnp,lwd=3,col='gray')
lines(t,fnp-3*errof,lty=2,lwd=3,col='gray')
lines(t,fnp+3*errof,lty=2,lwd=3,col='gray')


################################################################################

# 2. Real Data: ragweed pollen concentration
rm(list=ls(all=TRUE))
source('functions.R')

dados=read.table("ragweed.txt")
ind= order(dados[,3])
dado=dados[ind,]
n=nrow(dado)
t=dado[,3]# day.in.seas: non-parametric explanatory variable
y=sqrt(dado[,1])
year=dado[,2]

rain=dado[,5]
temp=dado[,4]
wind=dado[,6]

X=matrix(1,n,3)
X[,1]=rain
X[,2]=temp
X[,3]=wind
p=3

nk=1
for (i in 2:n){
    aux=t[i]-t[i-1]
    nki=i*(aux>0)
    nk=c(nk,nki)
}
t0=t[nk>0]#Knots

Nc=N(t,t0)
Rc=R(t0)
Qc=Q(t0)
K=Qc%*%solve(Rc)%*%t(Qc)
X1=X

# Estimating Normal and SN PLM 
tetaN=EMNspSIC_DF(y,X1,Nc,K)#theta=(beta,sigma^2,f,alpha)
miN=miNsp(tetaN,y,X1,Nc,K)
epN=sqrt(diag(solve(miN)))

teta=EMSNspSIC_DF(y,X1,Nc,K)# teta=(beta,sigma^2,lambda,f,alpha)
miSN=miSNsp(teta,y,X1,Nc,K)
epSN=sqrt(diag(solve(miSN))) 

logN=veroNsp(tetaN,y,X1,Nc,K)
logSNsp=veroSNsp(teta,y,X1,Nc,K)


# Criterions:
# Normal
alphan=as.numeric(tetaN[length(tetaN)])
sigma2n=as.numeric(tetaN[p+1])
Sf=solve(t(Nc)%*%Nc+alphan*sigma2n*K)%*%t(Nc)
glN=sum(diag(Nc%*%Sf))
pN= p+1+glN
AICN=-2*logN+2*pN
BICN=-2*logN+log(n)*pN

# SN
logSNsp=veroSNsp(teta,y,X1,Nc,K)
alphae=as.numeric(teta[length(teta)])
sigma2e=as.numeric(teta[p+1])
lambdae=as.numeric(teta[p+2])
Sf=solve(t(Nc)%*%Nc+alphae*sigma2e/(1+lambdae^2)*K)%*%t(Nc)
glSN=sum(diag(Nc%*%Sf))
pSN= p+2+glSN
AICSNsp=-2*logSNsp+2*pSN
BICSNsp=-2*logSNsp+log(n)*pSN

epfSN=epSN[(p+3):length(epSN)]# length(t0) x 1  
epfN=epN[(p+2):length(epN)]# length(t0) x 1
  

# Graphics of non-parametric component
# Normal
fnp=tetaN[(p+2):(length(tetaN)-1)]
bnp=tetaN[1:p]
yestNp=X%*%bnp
yestNnp=Nc%*%fnp
yestNd=y-X%*%bnp
bnup=fnp+2*epfN
bndo=fnp-2*epfN
plot(t,yestNd,type='p',xlab='Days in season',ylab='f(Days in season)',font.lab=2,ps=0.1)
lines(t,yestNnp,lwd=2,col='gray48')
lines(t0,bnup,lwd=2,lty=2,col='gray48')
lines(t0,bndo,lwd=2,lty=2,col='gray48')


# SN
p=ncol(X1)
bsnp=teta[1:p]
sigma2e=as.numeric(teta[p+1])
lambdae=as.numeric(teta[p+2])
vies=sqrt(2/pi)*sqrt(sigma2e)*lambdae/sqrt(1+lambdae^2)
fsnp=teta[(p+3):(length(teta)-1)]
yestSNd=y-(X1%*%bsnp+vies)
yestSNnp=Nc%*%fsnp# n x 1
bup=fsnp+2*epfSN
bdo=fsnp-2*epfSN
plot(t,yestSNd,type='p',xlab='Days in season',ylab='f(Days in season)',font.lab=2,ps=0.1)
lines(t,yestSNnp,lwd=2,,col='gray48')
lines(t0,bup,lwd=2,lty=2,col='gray48')
lines(t0,bdo,lwd=2,lty=2,col='gray48')


# Simulated envelopes
envelSNsp(y,X1,Nc,tetaN1,alpha=0.01)
envelSNsp(y,X1,Nc,teta,alpha=0.01)


####################################################################################
# Diagnostic analysis
V=3 # wind speed
aa=1:n

M0_casos=diag_Zhu(y,X1,teta,Nc,K,1,0,V) # y,X,theta,N,K,perturb,m0,V  , cases
corte=1/n+3*sd(M0_casos)
plot(aa,M0_casos,type='b',ylim=c(0,max(1.1*corte,max(M0_casos))),xlab='Index',ylab='M(0)',font.lab=2)
lines(aa,corte*rep(1,n),lty=2,lwd=2,col='grey50')
identify(aa,M0_casos,n=1)

M0_resp=diag_Zhu(y,X1,teta,Nc,K,2,0,V) # response
corte=1/n+3*sd(M0_resp)
plot(aa,M0_resp,type='b',ylim=c(0,max(1.1*corte,max(M0_resp))),xlab='Index',ylab='M(0)',font.lab=2)
lines(aa,corte*rep(1,n),lty=2,lwd=2,col='grey50')
identify(aa,M0_resp,n=2)

V=3# wind
M0_wind=diag_Zhu(y,X1,teta,Nc,K,3,0,V) 
corte=1/n+2*sd(M0_wind)
plot(aa,M0_wind,type='b',ylim=c(0,max(1.1*corte,max(M0_wind))),xlab='Index',ylab='M(0)',font.lab=2)
lines(aa,corte*rep(1,n),lty=2,lwd=2,col='grey50')

V=2# temp
M0_temp=diag_Zhu(y,X1,teta,Nc,K,3,0,V) 
corte=1/n+3*sd(M0_temp)
plot(aa,M0_temp,type='b',ylim=c(0,max(1.1*corte,max(M0_temp))),xlab='Index',ylab='M(0)',font.lab=2)
lines(aa,corte*rep(1,n),lty=2,lwd=2,col='grey50')
identify(aa,M0_temp,n=4)

Lever=Leverage(y,X1,teta,Nc,K)
corte=3*mean(Lever)
aa=1:n
plot(aa,Lever,type='b',ylim=c(0,max(1.1*corte,max(Lever))),xlab='Index',ylab='GL',font.lab=2)
lines(aa,corte*rep(1,n),lty=2,lwd=2,col='grey50')
identify(aa,Lever,n=3)

