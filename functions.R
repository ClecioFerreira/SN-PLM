require(moments)

auxpar<-function(y,X,Nc,K){
n=length(y)
p=ncol(X)
teptetam=matrix(0,p+2+p+3+n,1)
tetai=EMSNspSIC_DF(y,X,Nc,K) #[beta;sigma2;lambda;f;alpha]
#mi=miSNsp(tetai,y,X,Nc,K) # teta = [beta;sigma2;lambda;f;alpha]
#mi1=diag(solve(mi))
teptetam[1:(p+3+n)]=tetai
#teptetam[(p+3+n+1):(p+3+n+p+2)]=sqrt(mi1[1:3])
return(teptetam)
}

# Random samples
geraSN <-function(n,mu,sigma2,lambda){
delta=lambda/sqrt(1+lambda^2)
u1 = matrix(rnorm(n),n,1)
u2 = matrix(rnorm(n),n,1)
erro = sqrt(sigma2)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
y = mu+erro
return(y)
}

# Auxiliary functions

N <-function(t,t0){
n=length(t)
r=length(t0)
N=matrix(0,n,r)
for (i in 1:n){
    for (k in 1:r){
        if (t[i]==t0[k]) N[i,k]=1
    }
}
return(N)
}

Q <- function(t0){
r=length(t0)
h=rep(0,r-1)
for (i in 1:(r-1))  h[i]=t0[i+1]-t0[i]
q1=matrix(0,r,r-1)
for (i in 1:r){
   for (k in 2:(r-1)){
       if (abs(i-k)<2){
           q1[k-1,k]=1/h[k-1]
           q1[k,k]=-(1/h[k-1]+1/h[k])
           q1[k+1,k]=1/h[k]
       }
   }
}
Q=q1[1:r,2:(r-1)]
return(Q)
}

R <-function(t0){
r=length(t0)
h=rep(0,r-1)
for (i in 1:(r-1))  h[i]=t0[i+1]-t0[i]
r0=matrix(0,r-1,r-1)
for (i in 2:(r-2)){
    for (k in 1:(r-2)){
        if (abs(i-k)<2){
            r0[i,i]=1/3*(h[i-1]+h[i])
            r0[i,i+1]=h[i]/6
            r0[i+1,i]=h[i]/6
        }
    }
}
R=r0[2:(r-1),2:(r-1)]
return(R)
}


# EM algorithm

EMNspSIC_DF <-function(y,X,N,K){
alpha=optim(100,SICNspalpha,gr=NULL,y,X,N,K,method="L-BFGS-B",lower=0.001,upper=1000,control=list(maxit=50))# 0.001,1000
teta=EMNsp(y,X,N,K,alpha$par)
g=rbind(teta,alpha$par)
# theta=[beta,sigma2,f,alpha]
}

EMSNspSIC_DF <-function(y,X,N,K){
alpha=optim(10,SICSNspalpha,gr=NULL,y,X,N,K,method="L-BFGS-B",lower=0.001,upper=500,control=list(maxit=100))# 0.001,1000
teta=EMSNsp(y,X,N,K,alpha$par)
g=rbind(teta,alpha$par)
# theta=[beta,sigma2,f,alpha]
}


SICNspalpha<-function(alpha,y,X,N,K){
n=length(y)
p=ncol(X)
teta=EMNsp(y,X,N,K,alpha)#[beta,sigma2,f]
sigma2=as.numeric(teta[p+1])
logN=veroNsp(rbind(teta,alpha),y,X,N,K)
Sf=solve(t(N)%*%N+alpha*sigma2*K)%*%t(N)
SIC=-2*logN+(p+1+sum(diag(N%*%Sf)))*log(n)
}

SICSNspalpha<-function(alpha,y,X,N,K){
n=length(y)
p=ncol(X)
teta=EMSNsp(y,X,N,K,alpha)#[beta,sigma2,lambda,f]
sigma2=as.numeric(teta[p+1])
lambda=as.numeric(teta[p+2])
logSN=veroSNsp(rbind(teta,alpha),y,X,N,K)
Sf=solve(t(N)%*%N+alpha*sigma2/(1+lambda^2)*K)%*%t(N)
SIC=-2*logSN+(p+2+sum(diag(N%*%Sf)))*log(n)
}


EMNsp<-function(y,X,N,K,alpha){
n=length(y)
# valores iniciais
auxbeta=solve(t(X)%*%X)%*%t(X)
beta0=auxbeta%*%y
sigma2=as.numeric(t(y-X%*%beta0)%*%(y-X%*%beta0))/n
auxN=t(N)%*%(y-X%*%beta0)
f=solve(t(N)%*%N+alpha*sigma2*K)%*%auxN
theta0=rbind(beta0,sigma2)
criterio=norm(theta0)
cont=0
#print(alpha)
#print(sigma2)
while ((criterio > 1e-5)&&(cont<5000)){
    cont=cont+1
    beta1=auxbeta%*%(y-N%*%f)
    auxN=t(N)%*%(y-X%*%beta1)
    f=solve(t(N)%*%N+alpha*sigma2*K)%*%auxN
    sigma2=as.numeric(t(y-X%*%beta1-N%*%f)%*%(y-X%*%beta1-N%*%f))/n
    theta=rbind(beta1,sigma2)
    criterio=sqrt(sum((theta-theta0)^2))
    theta0=theta
}
g=rbind(theta,f)
return(g)
}

EMSNsp<-function(y,X,N,K,alpha){
n=length(y)
# valores iniciais
auxbeta=solve(t(X)%*%X)%*%t(X)
beta0=auxbeta%*%y
sigma2=as.numeric(t(y-X%*%beta0)%*%(y-X%*%beta0))/n
auxN=t(N)%*%(y-X%*%beta0)
f=solve(t(N)%*%N+alpha*sigma2*K)%*%auxN
lambda=as.numeric(skewness(y-X%*%beta0-N%*%f))
theta0=rbind(beta0,sigma2,lambda)
criterio=norm(theta0)
cont=0
#print(alpha)
#print(sigma2)
while ((criterio > 1e-5)&&(cont<5000)){
    cont=cont+1
    res=y-X%*%beta0-N%*%f
    sigma=sqrt(sigma2)
    eta=lambda*res
    aux=eta/sigma
    aux1=pmax(aux,-37)
    Wphi=dnorm(aux1)/pnorm(aux1)
    z=eta+sigma*Wphi
    z2=eta^2+sigma2+sigma*eta*Wphi
    beta0=auxbeta%*%(y-N%*%f-lambda/(1+lambda^2)*z)
    auxN=t(N)%*%(y-X%*%beta0-lambda/(1+lambda^2)*z)
    f=solve(t(N)%*%N+alpha*sigma2*K/(1+lambda^2))%*%auxN
    Qb=as.numeric(t(res)%*%res)
    sigma2=as.numeric((1+lambda^2)*Qb+sum(z2)-2*lambda*t(z)%*%res)/(2*n)
    lambda=as.numeric(t(z)%*%res)/Qb
    theta=rbind(beta0,sigma2,lambda)
    criterio=sqrt(sum((theta-theta0)^2))
    theta0=theta
}
g=rbind(theta,f)
return(g)
}

# Log-likelihoods
veroNsp<-function(theta,y,X,N,K){
# log-vero modelo normal semi-paramétrico
# theta=[beta;sigma2;f;alpha]
p=ncol(X)
beta0=theta[1:p]
sigma=as.numeric(sqrt(theta[p+1]))
f=theta[(p+2):(length(theta)-1)]
alpha=as.numeric(theta[length(theta)])
mu=X%*%beta0+N%*%f
logv=sum(log(dnorm(y,mu,sigma)))-alpha/2*t(f)%*%K%*%f
return(as.numeric(logv))
}

veroSNsp <-function(theta,y,X,N,K){
# theta=[beta;sigma2;lambda;f;alpha];
p=ncol(X)
beta0=theta[1:p]
sigma2=as.numeric(theta[p+1])
sigma=sqrt(sigma2)
lambda=as.numeric(theta[p+2])
f=theta[(p+3):(length(theta)-1)]
alpha=as.numeric(theta[length(theta)])
mu=X%*%beta0+N%*%f
vero=2*dnorm(y,mu,sigma)*pnorm(pmax(lambda*(y-mu)/sigma,-37))
g1=sum(log(vero))-alpha/2*t(f)%*%K%*%f
return(g1)
}

# Information Matrix
miSNsp <-function(teta,y,X,N,K){
# teta = [beta;sigma2;lambda;f;alpha];
n=length(y)
p=ncol(X)
q0=ncol(K)
beta0=teta[1:p]
sigma2=as.numeric(teta[p+1])
lambda=as.numeric(teta[p+2])
f=teta[(p+3):(length(teta)-1)]
alpha=as.numeric(teta[length(teta)])
res=y-X%*%beta0-N%*%f
sigma=sqrt(sigma2)
aux=lambda*res/sigma
aux1=pmax(aux,-37)
Wphi=dnorm(aux1)/pnorm(aux1)
######################### I2(beta,sigma2,lambda,f) ######################
Wphi1= as.numeric(-Wphi*(aux1+Wphi))
Qwbeta= as.numeric(t(res)%*%diag(Wphi1)%*%res)
mil=p+q0+2
I2=matrix(0,mil,mil)
I2[1:p,1:p]= -lambda^2/sigma2*t(X)%*%diag(Wphi1)%*%X
I2[1:p,p+1]= -lambda/(2*(sigma^3))*t(X)%*%Wphi-lambda^2/(2*(sigma2^2))*t(X)%*%diag(Wphi1)%*%res
I2[1:p,p+2]= 1/sigma*t(X)%*%Wphi+lambda/sigma2*t(X)%*%diag(Wphi1)%*%res
I2[p+1,p+1]= -3*lambda/(4*(sigma^5))*t(res)%*%Wphi-lambda^2/(4*(sigma2^3))*Qwbeta
I2[p+1,p+2]= 1/(2*(sigma^3))*t(res)%*%Wphi+lambda/(2*(sigma2^2))*Qwbeta
I2[p+2,p+2]= -Qwbeta/sigma2
I2[1:p,(p+3):mil]= -lambda^2/sigma2*t(X)%*%diag(Wphi1)%*%N
I2[(p+3):mil,1:p]= -lambda^2/sigma2*t(N)%*%diag(Wphi1)%*%X
I2[(p+3):mil,(p+3):mil]= -lambda^2/sigma2*t(N)%*%diag(Wphi1)%*%N
I2[(p+3):mil,p+1]= -lambda/(2*(sigma^3))*t(N)%*%Wphi-lambda^2/(2*(sigma2^2))*t(N)%*%diag(Wphi1)%*%res
I2[(p+3):mil,p+2]= 1/sigma*t(N)%*%Wphi+lambda/sigma2*t(N)%*%diag(Wphi1)%*%res
I2[p+2,p+1]= I2[p+1,p+2]
I2[p+1,1:p]= t(I2[1:p,p+1])
I2[p+2,1:p]= t(I2[1:p,p+2])
I2[(p+1),(p+3):mil]=t(I2[(p+3):mil,p+1])
I2[(p+2),(p+3):mil]=t(I2[(p+3):mil,p+2])
################### I1(beta,sigma2,lambda,f) ################## 
I1=matrix(0,mil,mil)  
I1[1:p,1:p]=1/sigma2*t(X)%*%X
I1[1:p,p+1]=1/(sigma2^2)*t(X)%*%res
I1[p+1,p+1]= -n/(2*(sigma2^2))+1/(sigma2^3)*t(res)%*%res
I1[p+1,1:p]=t(I1[1:p,p+1])
I1[(p+3):mil,(p+3):mil]=1/sigma2*t(N)%*%N+alpha*K
I1[1:p,(p+3):mil]= 1/sigma2*t(X)%*%N
I1[(p+3):mil,1:p]= 1/sigma2*t(N)%*%X
I1[(p+3):mil,(p+1)]=1/(sigma2^2)*t(N)%*%res
I1[(p+1),(p+3):mil]=t(I1[(p+3):mil,(p+1)])
g=I1+I2
return(g)
}

miNsp <-function(teta,y,X,N,K){
# teta = [beta;sigma2;f;alpha];
n=length(y)
p=ncol(X)
q0=ncol(K)
beta0=teta[1:p]
sigma2=as.numeric(teta[p+1])
f=teta[(p+2):(length(teta)-1)]
alpha=as.numeric(teta[length(teta)])
res=y-X%*%beta0-N%*%f
################### I1(beta,sigma2,f) ######################
mil=p+q0+1
I1=matrix(0,mil,mil)
I1[1:p,1:p]=1/sigma2*t(X)%*%X
I1[1:p,p+1]=1/(sigma2^2)*t(X)%*%res
I1[p+1,p+1]= -n/(2*(sigma2^2))+1/(sigma2^3)*t(res)%*%res
I1[p+1,1:p]=t(I1[1:p,p+1])
I1[(p+2):mil,(p+2):mil]=1/sigma2*t(N)%*%N+alpha*K
I1[1:p,(p+2):mil]= 1/sigma2*t(X)%*%N
I1[(p+2):mil,1:p]= 1/sigma2*t(N)%*%X
I1[(p+2):mil,(p+1)]=1/(sigma2^2)*t(N)%*%res
I1[(p+1),(p+2):mil]=t(I1[(p+2):mil,(p+1)])
return(I1)
}

###############################################################################################################
# SN linear regression model
snn.em<- function(y,X) {
 n<-length(y)
 if(missing(X)) {X<-matrix(1,n,1)}
 p<-ncol(X)
 lambda<-as.numeric(skewness(y))
 beta<-solve(t(X)%*%X)%*%t(X)%*%y
 sigma2<-as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-p))
 theta0<-rbind(beta,sigma2,lambda)
 criterio<-(sum(theta0^2))^0.5
 cont<-0
 while (criterio > 1e-6) {
    cont<-cont+1
    res<-y-X%*%beta
    sigma<-sqrt(sigma2)
    eta<-lambda*res
    aux<-eta/sigma
    aux1<-pmax(aux,-37)
    Wphi<-dnorm(aux1)/pnorm(aux1)
    t1<-eta+sigma*Wphi
    t2<-eta^2+sigma2+sigma*eta*Wphi
    d<-res^2/sigma2
    beta<-solve(t(X)%*%(diag(1,n)+lambda^2*diag(1,n))%*%X)%*%t(X)%*%(y-lambda*(t1-lambda*y))
    Qb<-t(res)%*%res
    sigma2<-as.numeric((Qb+sum(t2)-2*lambda*t(t1)%*%res+lambda^2*Qb)/(2*n))
    lambda<-as.numeric(t(t1)%*%res/Qb)
    theta<-rbind(beta,sigma2,lambda)
    dif<-theta-theta0
    criterio<-(sum(dif^2))^0.5
    theta0<-theta
 }
  return(theta)
 }

 sn.logL <- function (theta,y,X) 
{
  #% Theta=[beta,sigma2,lambda]
  p<-ncol(X)
  beta<-theta[1:p]
  sigma2<-theta[p+1]
  lambda<-theta[p+2]
  z=(y-X%*%beta)/sqrt(sigma2)
  d<-z^2
  aux<-signif(lambda*z,digits=7)
  aux1<-signif(pmax(aux,-36),digits=7)
  cte<-2/sqrt(2*pi*sigma2)
  vero<-cte*exp(-0.5*d)*pnorm(aux1)
  return(sum(log(vero)))
}

##############################################################################################################
# Simulated envelopes

envelSNsp <- function(y,X,N,theta,alpha=0.05){
  n=length(y)
  p=ncol(X)
  beta0 = theta[1:p]
  sigma2=as.numeric(theta[p+1])
  f=theta[(p+3):(length(theta)-1)]
  mu=X%*%beta0+N%*%f
  d2=as.numeric((y-mu)^2/sigma2)
  d2s=sort(d2)
  d2s=t(d2s)
  xq2 <- qchisq(ppoints(n), 1)
  replic=200
  Xsim<-matrix(0,replic,n)
  for(i in 1:replic) Xsim[i,]<-rchisq(n, 1)
  Xsim2<-apply(Xsim,1,sort)
  d21<-rep(0,n)
  d22<-rep(0,n)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,],alpha/2)
    d22[i]  <- quantile(Xsim2[i,],1-alpha/2)}
  d2med <-apply(Xsim2,1,mean)
  fy <- range(d2s,d21,d22)
  
aa<-sort(d2,index.return=TRUE)
tsi2<-aa$x
ordem=aa$ix

aux<-0
for(i in 1:n){
if ((tsi2[i]<d21[i])|(tsi2[i]>d22[i])) aux<-c(aux,i) # |=ou
}
#print(aux)
ll=length(aux)
aux2=aux[2:ll]
if (ll >= 1) posi=ordem[aux2]
  
  plot(xq2,d2s,xlab = expression(bold(paste("Theoretical ",chi[1]^2, " quantile"))),
       ylab="Sample value",pch=20,ylim=fy,font.lab=2)
  par(new=T)
  plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
  par(new=T)
  plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="",lty="dashed")
  par(new=T)
  plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
  print(posi)
}


envelSNp <- function(y,X,theta,alpha=0.05){
  n=length(y)
  p=ncol(X)
  mu=X%*%theta[1:p]
  d2=as.numeric((y-mu)^2/theta[p+1])
  d2s=sort(d2)
  d2s=t(d2s)
  xq2 <- qchisq(ppoints(n), 1)
  replic=200
  Xsim<-matrix(0,replic,n)
  for(i in 1:replic) Xsim[i,]<-rchisq(n, 1)
  Xsim2<-apply(Xsim,1,sort)
  d21<-rep(0,n)
  d22<-rep(0,n)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,],alpha/2)
    d22[i]  <- quantile(Xsim2[i,],1-alpha/2)}
  d2med <-apply(Xsim2,1,mean)
  fy <- range(d2s,d21,d22)
  plot(xq2,d2s,xlab = expression(bold(paste("Theoretical ",chi[1]^2, " quantiles"))),
       ylab="Sample values and simulated envelope",pch=20,ylim=fy,font.lab=2)
  par(new=T)
  plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
  par(new=T)
  plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="",lty="dashed")
  par(new=T)
  plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
}

envelNsp<-function(y,X,N,K,theta){
#theta=(beta,sigma^2,f,alpha)  
n <- nrow(X)
p <- ncol(X)
beta0 = theta[1:p]
sigma2=as.numeric(theta[p+1])
f=theta[(p+2):(length(theta)-1)]
alpha=as.numeric(theta[length(theta)])
q0=ncol(N)
pq=p+q0;
MXN=matrix(0,pq,pq)
MXN[1:p,1:p]=t(X)%*%X
MXN[1:p,(p+1):pq]=t(X)%*%N
MXN[(p+1):pq,1:p]=t(N)%*%X
MXN[(p+1):pq,(p+1):pq]=t(N)%*%N+alpha*K
Minv=solve(MXN);
aux=matrix(0,n,pq)
aux[,1:p]=X
aux[,(p+1):pq]=N
H=aux%*%Minv%*%t(aux)
varR=sigma2*(diag(n)-H)%*%t(diag(n)-H)
tsi=as.numeric((diag(n)-H)%*%y)/sqrt(diag(varR))
#print(tsi)
#
aa<-sort(tsi,index.return=TRUE)
tsi2<-aa$x
ordem=aa$ix
#
ident <- diag(n)
epsilon <- matrix(0,n,100)
e <- matrix(0,n,100)
e1 <- numeric(n)
e2 <- numeric(n)
#
for(i in 1:100){
     epsilon[,i] <- rnorm(n,0,1)
     e[,i] <- (ident - H)%*%epsilon[,i]
     u <- diag(ident - H)
     e[,i] <- e[,i]/sqrt(u)
     e[,i] <- sort(e[,i]) }
#
for(i in 1:n){
     eo <- sort(e[i,])
     e1[i] <- (eo[2]+eo[3])/2
     e2[i] <- (eo[97]+eo[98])/2 }
#
med <- apply(e,1,mean)
faixa <- range(tsi,e1,e2)
#
aux<-0
for(i in 1:n){
if ((tsi2[i]<e1[i])|(tsi2[i]>e2[i])) aux<-c(aux,i) # |=ou
}
#print(aux)
ll=length(aux)
aux2=aux[2:ll]
if (ll >= 1) posi=ordem[aux2]
#
par(pty="s")
qqnorm(tsi,xlab="Quantiles of the standard normal",
ylab="Studentized residuals", ylim=faixa, pch=16,font.lab=2)
par(new=T)
qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1)
par(new=T)
qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1)
par(new=T)
qqnorm(med,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=2)
#--------------------------------------------------------------#
envel2<-posi
print(posi)

}


envelNsp_gam<-function(y,X,N,t,alpha){
p=ncol(X)
q0=ncol(N)
pq=p+q0;
MXN=matrix(0,pq,pq)
MXN[1:p,1:p]=t(X)%*%X
MXN[1:p,(p+1):pq]=t(X)%*%N
MXN[(p+1):pq,1:p]=t(N)%*%X
MXN[(p+1):pq,(p+1):pq]=t(N)%*%N+alpha*K
Minv=solve(MXN)
aux=matrix(0,n,pq)
aux[,1:p]=X
aux[,(p+1):pq]=N
H=aux%*%Minv%*%t(aux) 
n=length(fit.model$y)
mu=fit.model$linear.predictors
phi=fit.model$sig2
varr=as.numeric(diag(diag(n)-H))*phi
ri=fit.model$residuals
tdf=ri/sqrt(varr)
e=matrix(0,n,100)
resp=NULL
for(i in 1:100){
     resp<-rnorm(n,0,1)
     resp<-mu+sqrt(phi)*resp
     fit<-gam(resp~X[,1]+X[,2]+s(t,bs="cs", fx=FALSE), method="GCV.Cp", family=gaussian)
     phi=fit$sig2
     varr=as.numeric(diag(diag(n)-H))*phi
     ri<-fit$residuals
     td<-ri/sqrt(varr)
     e[,i]=sort(td)
}
e1 <- numeric(n)
e2 <- numeric(n)
e8 <- numeric(n)
for(i in 1:n){
     eo <- sort(e[i,])	
     e1[i] <- eo[ceiling(100*0.025)]
     e8[i] <- eo[ceiling(100*0.50)]
     e2[i] <- eo[ceiling(100*0.975)]
}
faixa <- range(tdf,e1,e2)
par(pty="s")
qqnorm(tdf,xlab="Percentiles of N(0,1)",
ylab="Standard Residuals", ylim=faixa, pch=16,font.lab=2)
par(new=TRUE)
qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1)
par(new=TRUE)
qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1)
par(new=TRUE)
qqnorm(e8,axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2)   
}

#############################################################################################################
# Local influence
diag_Zhu <- function(y,X,theta,N,K,perturb,m0,V) {
n=length(y)
p=ncol(X)
######## Delta's 
Q2=Hessian(y,X,theta,N,K)
if (perturb==1) delta=cases(y,X,theta,N,K) # casos
else if (perturb==2) delta=response(y,X,theta,N,K) # resposta
else if (perturb==3) delta=explic(y,X,theta,N,K,V) # explicativa
else delta=0
cd=2*t(delta)%*%solve(-Q2)%*%delta # size=pert x pert
cd1=(cd+t(cd))/2
AVE=eigen(cd1)
H=AVE$vectors
y1=AVE$values #dah os autovalores, VETOR (ordem decrescente): n x 1
r=sum(y1>0)
# Autovalores em ordem crescente autovetores na mesma ordem por colunas.
# H(i,:)*H(i,:)'=1, para todo i
ind1=y1>1e-10
ind2=y1>=m0
ind3=ind1*ind2#[ind1 ind2 ind3 y1]
r2=which(ind3==1) # posiçao
lambdapad=y1[r2]/r
Mm0=matrix(0,n,1)
H2=H*H
cont=0
for (i in 1:min(r2)){
    cont=cont+1
    aux=lambdapad[cont]*H2[,i]
    Mm0=Mm0+aux
}
return(Mm0)
}

cases <-function(y,X,theta,N,K)  {
n=length(y)
p=ncol(X)
beta1=theta[1:p]
sigma2=as.numeric(theta[p+1])
sigma=sqrt(sigma2)
lambda=as.numeric(theta[p+2])
f=theta[(p+3):(length(theta)-1)]
alpha=as.numeric(theta[length(theta)])
res=y-X%*%beta1-N%*%f
Sbetaf=as.numeric(t(res)%*%res)
aux    <- lambda*res/sigma
aux1   <- pmax(aux,-37)
Wphi   <- dnorm(aux1)/pnorm(aux1)
#d      <- res^2/sigma2
z     <- lambda*res+sigma*Wphi
z2     <- lambda^2*res^2+sigma2+lambda*sigma*res*Wphi
pend= p+2+length(f)
Delta=matrix(0,pend,n)
Delta[1:p,]=1/sigma2*t(X)%*%diag(-lambda*as.vector(z)+(1+lambda^2)*as.vector(res))
Delta[(p+1),]=-1/sigma2*rep(1,n)+1/(2*(sigma2^2))*((1+lambda^2)*t(res)%*%diag(as.vector(res))+t(z2)-2*lambda*t(res)%*%diag(as.vector(z)))
Delta[(p+2),]=1/sigma2*t(res)%*%diag(as.vector(z)-lambda*as.vector(res))
Delta[(p+3):pend,]=1/sigma2*t(N)%*%diag(-lambda*as.vector(z)+(1+lambda^2)*as.vector(res))
#print(Delta)
return(Delta)
}

response <-function(y,X,theta,N,K) {
n=length(y)
p=ncol(X)
beta1=theta[1:p]
sigma2=as.numeric(theta[p+1])
sigma=sqrt(sigma2)
lambda=as.numeric(theta[p+2])
f=theta[(p+3):(length(theta)-1)]
alpha=as.numeric(theta[length(theta)])
res=y-X%*%beta1-N%*%f
Sbetaf=as.numeric(t(res)%*%res)
aux    <- lambda*res/sigma
aux1   <- pmax(aux,-37)
Wphi   <- dnorm(aux1)/pnorm(aux1)
#d      <- res^2/sigma2
z     <- lambda*res+sigma*Wphi
z2     <- lambda^2*res^2+sigma2+lambda*sigma*res*Wphi
pend= p+2+length(f)
Sy=sd(y)
Delta=matrix(0,pend,n)
Delta[1:p,]=(1+lambda^2)*t(X)
Delta[(p+1),]=1/sigma2*((1+lambda^2)*t(res)-lambda*t(z))
Delta[(p+2),]=t(z)-2*lambda*t(res)
Delta[(p+3):pend,]= (1+lambda^2)*t(N)
return(Sy/sigma2*Delta)
}

explic <-function(y,X,theta,N,K,V) {
n=length(y)
p=ncol(X)
beta1=theta[1:p]
sigma2=as.numeric(theta[p+1])
sigma=sqrt(sigma2)
lambda=as.numeric(theta[p+2])
f=theta[(p+3):(length(theta)-1)]
alpha=as.numeric(theta[length(theta)])
res=y-X%*%beta1-N%*%f
Sbetaf=as.numeric(t(res)%*%res)
aux    <- lambda*res/sigma
aux1   <- pmax(aux,-37)
Wphi   <- dnorm(aux1)/pnorm(aux1)
#d      <- res^2/sigma2
z     <- lambda*res+sigma*Wphi
z2     <- lambda^2*res^2+sigma2+lambda*sigma*res*Wphi
pend= p+2+length(f)
n1=matrix(0,p,1)
n1[V]=1
betat=as.numeric(beta1[V])
St=sd(X[,V])
Delta=matrix(0,pend,n)
Delta[1:p,]=(1+lambda^2)*(n1%*%t(res)-betat*t(X))-lambda*n1%*%t(z)
Delta[(p+1),]=betat/sigma2*(lambda*t(z)-(1+lambda^2)*t(res))
Delta[(p+2),]=betat*(-t(z)+2*lambda*t(res))
Delta[(p+3):pend,]= -betat*(1+lambda^2)*t(N)
return(St/sigma2*Delta)
}

Leverage <-function(y,X,theta,N,K) {
n=length(y)
p=ncol(X)
beta1=theta[1:p]
sigma2=as.numeric(theta[p+1])
sigma=sqrt(sigma2)
lambda=as.numeric(theta[p+2])
delta=lambda/sqrt(1+lambda^2)
f=theta[(p+3):(length(theta)-1)]
alpha=as.numeric(theta[length(theta)])
res=y-X%*%beta1-N%*%f
aux    <- lambda*res/sigma
aux1   <- pmax(aux,-37)
Wphi   <- dnorm(aux1)/pnorm(aux1)
#d      <- res^2/sigma2
z     <- lambda*res+sigma*Wphi
z2     <- lambda^2*res^2+sigma2+lambda*sigma*res*Wphi
pend=p+2+length(f)
Dmu=matrix(0,n,pend)
Dmu[,1:p]=X
Dmu[,p+1]=0.5*sqrt(2/pi)*delta/sqrt(sigma2)
Dmu[,p+2]=sqrt(2/pi)*sigma/((1+lambda^2)^1.5)
Dmu[,(p+3):pend]=N
Delta=matrix(0,pend,n)
Delta[1:p,]=(1+lambda^2)*t(X)
Delta[(p+1),]=1/sigma2*((1+lambda^2)*t(res)-lambda*t(z))
Delta[(p+2),]=t(z)-2*lambda*t(res)
Delta[(p+3):pend,]= (1+lambda^2)*t(N)
Qtetay=1/sigma2*Delta
Qinv=solve(-Hessian(y,X,theta,N,K))
GL=diag(Dmu%*%Qinv%*%Qtetay)
return(GL)
}

Hessian<- function(y,X,theta,N,K){
# Matriz Hessiana: der^2 Q(theta,theta)/dtheta^T dtheta
# teta: (beta,sigma^2,eta,tau^2)
n=length(y)
p=ncol(X)
beta1=theta[1:p]
sigma2=as.numeric(theta[p+1])
sigma=sqrt(sigma2)
lambda=as.numeric(theta[p+2])
f=theta[(p+3):(length(theta)-1)]
alpha=as.numeric(theta[length(theta)])
res=y-X%*%beta1-N%*%f
Sbetaf=as.numeric(t(res)%*%res)
aux    <- lambda*res/sigma
aux1   <- pmax(aux,-37)
Wphi   <- dnorm(aux1)/pnorm(aux1)
#d      <- res^2/sigma2
z     <- lambda*res+sigma*Wphi
z2     <- lambda^2*res^2+sigma2+lambda*sigma*res*Wphi
pend= p+2+length(f)
H=matrix(0,pend,pend)
H[1:p,1:p]=-(1+lambda^2)/sigma2*t(X)%*%X
H[1:p,(p+1)]=-1/(sigma2^2)*t(X)%*%((1+lambda^2)*res-lambda*z)
H[1:p,(p+2)]= 1/sigma2*t(X)%*%(2*lambda*res-z)
H[1:p,(p+3):pend]= -(1+lambda^2)/sigma2*t(X)%*%N
H[(p+1),1:p]=t(H[1:p,(p+1)])
H[(p+1),(p+1)]=n/(sigma2^2)-1/(sigma2^3)*(sum(z2)-2*lambda*t(res)%*%z+(1+lambda^2)*Sbetaf)
H[(p+1),(p+2)]=1/(sigma2^2)*(lambda*Sbetaf-t(res)%*%z)
H[(p+1),(p+3):pend] = -1/(sigma2^2)*t(t(N)%*%((1+lambda^2)*res-lambda*z))
H[(p+2),1:p]=t(H[1:p,(p+2)])
H[(p+2),(p+1)]= H[(p+1),(p+2)]
H[(p+2),(p+2)]=-1/sigma2*Sbetaf
H[(p+2),(p+3):pend]= 1/sigma2*t(t(N)%*%(2*lambda*res-z))
H[(p+3):pend,1:p]=t(H[1:p,(p+3):pend])
H[(p+3):pend,(p+1)] = t(H[(p+1),(p+3):pend])
H[(p+3):pend,(p+2)]=t(H[(p+2),(p+3):pend])
H[(p+3):pend,(p+3):pend]= -(1+lambda^2)/sigma2*t(N)%*%N-alpha*K
return(H)
}
