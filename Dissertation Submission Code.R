###########################################################
###         MODULARIZATION IN BAYESIAN STATISTICS       ###
### a study on the cut model and semi-modular inference ###
###########################################################

# Code for the dissertation written by Hubert Pawlusinski 2011320 (MMORSE)
# University of Warwick, 2024













library(MASS)
if (!require(mvtnorm)) install.packages('mvtnorm')
library(mvtnorm)
library(rlist)


##################
### CONVENTION ###
##################

#For all prior and likelihood functions, argument for theta must be th and phi for phi.
#However! Contrary to the dissertation, theta represents the parameter of the first module
#         and phi represents a parameter of the second module, so parameters are swapped.
#         This will be changed eventually for the purpose of producing plots, but the sampler
#         and a lot of operations for the sampler use the other naming convention.
#For likelihoods Z and Y must be the names of the arguments, same as in the dissertation.

####################
### MH PROPOSALS ###
####################

#for Metropolis-Hastings proposals, we set the following:

MVN_S_prop=function(d,mean=NULL,scale){
  if (is.null(mean)){
    mean=rep(0,times=d)
  }
  sigma=diag(x=scale, nrow=d, ncol=d)
  prop=mvrnorm(n=1,mu=mean,Sigma=sigma)
  p=dmvnorm(x=prop,mean=mean,sigma=sigma)
  
  return(list(prop,p))
}

#####################
### FULL BAYESIAN ###
#####################
#VERY CONFUSINGLY I NAMED 1ST MODULE PARAMETER THETA AND 2ND MODULE PARAMETER PHI

#edit 1: we changed all comparisons onto logs - easier to compute and reduced risk of underflow
Bayes_SRW_sampler = function(N=1000, log_pr_th, log_pr_phi, l_Z, l_Y,
                             th0=NULL,phi0=NULL,prop,
                             phi_names,th_names,
                             Z,Y,scale1,scale2){ #no means as we perform SRW MH
  
  d1=length(th_names)
  d2=length(phi_names)
  thetas=data.frame(rep(list(NA),d1))
  colnames(thetas)=th_names
  phis=data.frame(rep(list(NA),d2))
  colnames(phis)=phi_names
  
  acc=c()
  
  if(is.null(th0)){
    th0=rep(0,d1)
  }
  thetas[1,]=th0
  if(is.null(phi0)){
    phi0=rep(0,d2)
  }
  phis[1,]=phi0
  
  th=th0
  phi=phi0
  #first sampler
  U=runif(N)
  for (i in 1:N){
    acc[i]=0
    prop1=prop(d=d1,scale=scale1)
    prop2=prop(d=d2,scale=scale2)
    e1=prop1[[1]]
    e2=prop2[[1]] #AAAAAA
    # print("e")
    # print(e1)
    # pq=prop[[2]] not needed for SRW Metropolis
    log_numer=log_pr_th(th+e1)+log_pr_phi(phi=phi+e2,th=th+e1)+l_Z(Z=Z,th=th+e1)+l_Y(Y=Y,phi=phi+e2,th=th+e1)
    log_denom=log_pr_th(th)+log_pr_phi(phi=phi,th=th)+l_Z(Z=Z,th=th)+l_Y(Y=Y,phi=phi,th=th)
    log_accrate=log_numer-log_denom
    # print("numer")
    # print(numer)
    # print("denom")
    # print(denom)
    # print("accrate")
    # print(accrate)
    if(U[i]<exp(log_accrate)){
      th=th+e1
      phi=phi+e2
      acc[i]=1
    }
    thetas[i+1,]=th
    phis[i+1,]=phi
  }
  
  return(list(thetas,phis,acc))
}

##########################
### CUT NESTED SAMPLER ###
##########################


#namings will follow papers from which are inspired, e.g. Carmona-Nicholls - CN sampler
CN_cut_sampler = function(N1=1000,N2=20, log_pr_th, log_pr_phi, l_Z, l_Y,
                          th0=NULL,phi0=NULL,prop1,prop2,
                          phi_names,th_names,
                          Z,Y,scale1,scale2){ #no means as we perform SRW MH
  
  d1=length(th_names)
  d2=length(phi_names)
  thetas=data.frame(rep(list(NA),d1))
  colnames(thetas)=th_names
  phis=data.frame(rep(list(NA),d2))
  colnames(phis)=phi_names
  
  acc1=c()
  acc2=c()
  
  if(is.null(th0)){
    th0=rep(0,d1)
  }
  thetas[1,]=th0
  if(is.null(phi0)){
    phi0=rep(0,d2)
  }
  phis[1,]=0
  
  th=th0
  #first sampler
  U=runif(N1)
  for (i in 1:N1){
    acc1[i]=0
    prop=prop1(d=d1,scale=scale1)
    e=prop[[1]]
    # pq=prop[[2]] not needed for SRW Metropolis
    log_accrate=log_pr_th(th+e)+l_Z(Z=Z,th=th+e)-log_pr_th(th)-l_Z(Z=Z,th=th)
    if(U[i]<exp(log_accrate)){
      th=th+e
      acc1[i]=1
    }
    thetas[i+1,]=th
  }
  
  #second sampler - repeated for each entry of the first one
  for (i in 1:N1){
    U=runif(N2)
    subacc=c()
    theta=thetas[i,]
    phi=phi0
    for (j in 1:N2){
      subacc[j]=0
      prop=prop2(d=d2,scale=scale2)
      e=prop[[1]]
      # pq=prop[[2]] not needed for SRW Metropolis
      log_accrate=log_pr_phi(phi=phi+e,th=theta)+l_Y(Y=Y,th=theta,phi=phi+e)-log_pr_phi(phi=phi,th=theta)-l_Y(Y=Y,th=theta,phi=phi)
      if(U[j]<exp(log_accrate)){
        phi=phi+e
        subacc[j]=1
      }
    }
    phis[i+1,]=phi
    acc2[i]=mean(subacc)
  }
  return(list(thetas,phis,acc1,acc2))
}

##########################
### SMI NESTED SAMPLER ###
##########################


CN_smi_sampler = function(N1=1000,N2=20, log_pr_th, log_pr_phi, l_Z, l_Y,
                          th0=NULL,phi0=NULL,prop1,prop2,
                          phi_names,th_names,
                          Z,Y,scale1,scale2, 
                          eta){ #no means as we perform SRW MH
  
  #setup of samples from data frames
  d1=length(th_names)
  d2=length(phi_names)
  thetas=data.frame(rep(list(NA),d1))
  colnames(thetas)=th_names
  phis=data.frame(rep(list(NA),d2))
  colnames(phis)=phi_names
  
  acc1=c()
  acc2=c()
  esjd=c()
  
  #initialization of the chain
  if(is.null(th0)){
    th0=rep(0,d1)
  }
  thetas[1,]=th0
  if(is.null(phi0)){
    phi0=rep(0,d2)
  }
  phis[1,]=phi0
  
  
  th=th0
  auphi=phi0
  #outer sampler
  U=runif(N1)
  for (i in 1:N1){
    acc1[i]=0
    prop=prop1(d=d1,scale=scale1)
    e1=prop[[1]]
    prop=prop1(d=d2,scale=scale2)
    e2=prop[[1]]
    # pq=prop[[2]] probability of the proposal not needed for SRW Metropolis
    #acceptance probability is computed on the logs for computational efficiency
    log_accrate_num=l_Z(Z,th=th+e1)+eta*l_Y(Y=Y,th=th+e1,phi=auphi+e2)+log_pr_th(th+e1)+log_pr_phi(phi=auphi+e2,th=th+e1)
    log_accrate_denom=l_Z(Z,th=th)+eta*l_Y(Y=Y,th=th,phi=auphi)+log_pr_th(th)+log_pr_phi(phi=auphi,th=th)
    if(U[i]<exp(log_accrate_num-log_accrate_denom)){
      th=th+e1
      auphi=auphi+e2
      acc1[i]=1
    }
    thetas[i+1,]=th
  }
  
  #inner sampler - repeated for each entry of the first one
  for (i in 1:N1){
    U=runif(N2)
    subacc=c()
    theta=thetas[i,]
    phi=phi0
    for (j in 1:N2){
      subacc[j]=0
      prop=prop2(d=d2,scale=scale2)
      e=prop[[1]]
      # pq=prop[[2]] proposal density not needed for SRW Metropolis
      log_accrate=log_pr_phi(phi=phi+e,th=theta)+l_Y(Y=Y,th=theta,phi=phi+e)-log_pr_phi(phi=phi,th=theta)-l_Y(Y=Y,th=theta,phi=phi)
      if(U[j]<exp(log_accrate)){
        phi=phi+e
        subacc[j]=1
      }
    }
    phis[i+1,]=phi
    sjd=sum((phis[i,]-phis[i+1,])^2)+sum((thetas[i,]-thetas[i+1,])^2)
    esjd[i]=sjd
    acc2[i]=mean(subacc)
  }
  return(list(thetas,phis,acc1,acc2,esjd))
}

####################
### ENERGY SCORE ###
####################
#we take sample of predictions of a vector Y and compare it with observation y from Y
#here we assume sample is iid from Y, which is important for computation of e2
#b is the value of beta for energy score
energy_1D=function(sample,y,b){
  l=length(sample)
  l_y=length(y)
  #compare sample with value y
  e1=0
  for (j in 1:l_y){
    for(i in 1:l){
      e1=e1+(abs(sample[i]-y[j]))^b
    }
  }
  e1=e1/l
  #divide sample into two
  d=floor(l/2)
  sample1=sample[1:d]
  sample2=sample[(d+1):(2*d)]
  e2=0
  #compare two 'independent' samples of X
  for (i in 1:d){
    e2=e2+(abs(sample1[i]-sample2[i]))^b
  }
  e2=e2*l_y/d
  return(2*e1-e2)
}


# We need to set up necessary likelihoods and priors
# log prior in first module
# log prior in 2nd module
# log likelihood of first module
# log likelihood of 2nd module


###################################
### BIASED DATA EXAMPLE (CH. 5) ###
###################################

### simulated data ####
set.seed(2011320)
n=25
m=50

Z=rnorm(n=n,mean=0,sd=2)
Y=rnorm(n=m,mean=1,sd=1)


#####

### probabilities for inference ####
#1st log prior
log_prior_th=function(th){
  return(log(dnorm(th,mean=0,sd=1)))
}
#2nd log prior
log_prior_phi=function(phi,th){
  return(log(dnorm(phi,mean=0,sd=0.1)))
}
#1st log likelihood
L_Z=function(Z,th){
  l=length(Z)
  p=1
  for (i in 1:l){
    p=p*dnorm(x=Z[i],mean=th,sd=2)
  }
  return(p)
}
#2nd log likelihood
l_Z=function(Z,th,sd=2){
  n=length(Z)
  p=(-n/2)*log(2*pi)-n*log(sd)-(2*sd^2)^(-1)*(sum((Z-th)*(Z-th)))
  return(p)
}
#likelihoods (not used as we operate on log likelihoods)
L_Y=function(Y,th,phi){
  l=length(Y)
  p=1
  for (i in 1:l){
    p=p*dnorm(x=Y[i],mean=th+phi,sd=1)
  }
  return(p)
}

l_Y=function(Y,th,phi,sd=1){
  n=length(Y)
  p=(-n/2)*log(2*pi)-n*log(sd)-(2*sd^2)^(-1)*(sum((Y-th-phi)*(Y-th-phi)))
  return(p)
}


#### Generating numbers for scoring ###
gen_Z=function(th){
  return(rnorm(n=1,mean=th,sd=2))
}

gen_Y=function(th,phi){
  return(rnorm(n=1,mean=th+phi,sd=1))
}

#####


### Choosing scale (Sec. 5.2.1.) #####
eta=0
scs=c(0.0025,0.01,0.04,0.16,0.25,1,4,25)
jumps=c()
for (k in 1:length(scs)){
  sc=scs[k]
  sample_smi=CN_smi_sampler(N1=1000,N2=30, log_pr_th=log_prior_th, log_pr_phi=log_prior_phi, l_Z=l_Z, l_Y=l_Y,
                            th0=NULL,phi0=NULL,prop1=MVN_S_prop,prop2=MVN_S_prop,
                            phi_names=c("phi"),th_names=c("theta"),
                            Z=Z,Y=Y,scale1=sc,scale2=0.1, eta=eta)
  jumps[k]=mean(sample_smi[[5]][-(1:100)])
}
jumps
max(jumps)
scs[which(jumps==max(jumps))] 
#0.16 optimal for standard training set
#0.25 optimal for double, or eta0.2
#x8 ranges from 0.04 to 0.25 so we pick 0.16
#1 for eta=0, 0.1 for small datasets
#so 1 for eta 0,0.1 for 0.4, 1, 2, 4 and 0.16 for everything else.

#####

### Adjusting inner sampler length (Sec. 5.2.2.) #####

#sampling from (theta|Y,phi) for a range of phis

phis=seq(-0.4,1,0.2)
ms=c(5,10,15,20,30,50,100)
N=1000
scale=0.4

samp=array(0,dim=c(length(phis),length(ms),N),dimnames=list(as.character(phis),as.character(ms)))

for (phi in phis){
  acc=c()
  for (i in 1:N){
    X=c()
    x=0
    U=runif(max(ms))
    for (j in 1:max(ms)){
      acc[(i-1)*ms+j]=0
      e=MVN_S_prop(d=1,scale=scale)[[1]]
      log_numer=log_prior_th(x+e)+l_Y(Y=Y,phi=phi,th=x+e)
      log_denom=log_prior_th(x)+l_Y(Y=Y,phi=phi,th=x)
      log_accrate=log_numer-log_denom
      if(log_accrate>log(U[j])){
        x=x+e
        acc[(i-1)*ms+j]=1
      }
      X[j]=x
    }
    for (m in ms){
      samp[as.character(phi),as.character(m),i]=X[m]
    }
  }
  print(phi)
  print(mean(acc))
}


#scale 0.1 yields very nice accuracy rates


#plots

phi=as.character(phis[1])

cols=rainbow(length(ms))


for (i in 1:length(ms)){
  m=ms[i]
  x=density(samp[phi,as.character(m),])
  if(i==1){
    plot(x,type="l",col=cols[i],ylim=c(0,3),main=paste("phi =",phi), xlab="theta")
  }else{
    points(x,type="l",col=cols[i])
  }
}
legend(-0.2,3,lty=1,legend=ms,col=cols, title="MCMC step")

#20 is the safe inner chain length as it already almost converges to what it should be
#30 for case when 4x as much data (not a huge difference but just to be safe)
#15 is good for scale 0.4 instead of 0.1 for smallest set


#####

######## Qualitative samplers checks (Sec. 5.3.) #######

### Full Bayesian ####
set.seed(2011320)
sample_FB=Bayes_SRW_sampler(N=10000,log_pr_th=log_prior_th,log_pr_phi=log_prior_phi,
                            l_Z=l_Z,l_Y=l_Y,
                            th0=NULL, phi0=NULL, prop=MVN_S_prop, 
                            phi_names = c("phi"), th_names = c("theta"),
                            Z=Z, Y=Y, scale1=0.025, scale2=0.1)
par(mfrow=c(5,2))
thetas_FB=sample_FB[[1]]
phis_FB=sample_FB[[2]]
hist(thetas_FB[,1],breaks=seq(from=-2,to=2,by=0.05),main="True mean for full Bayesian", xlab="phi")
abline(v=0,lty="dashed",col="red")
hist(phis_FB[,1],breaks=seq(from=-2,to=2,by=0.05),main="Bias for full Bayesian", xlab="theta")
abline(v=1,lty="dashed",col="red")
mean(sample_FB[[3]])


#####


### Cut model using CN nested MCMC ####
set.seed(2011320)
sample_cut=CN_cut_sampler(N1=3000,N2=20, log_pr_th=log_prior_th, log_pr_phi=log_prior_phi, l_Z=l_Z, l_Y=l_Y,
                          th0=NULL,phi0=NULL,prop1=MVN_S_prop,prop2=MVN_S_prop,
                          phi_names=c("phi"),th_names=c("theta"),
                          Z=Z,Y=Y,scale1=0.25,scale2=0.1)

phis_cut=sample_cut[[2]]
thetas_cut=sample_cut[[1]]
hist(thetas_cut[,1],breaks=seq(from=-2,to=2,by=0.05),main="True mean for cut posterior", xlab="phi")
abline(v=0,lty="dashed",col="red")
hist(phis_cut[,1],breaks=seq(from=-2,to=2,by=0.05),main="Bias for cut posterior", xlab="theta")
abline(v=1,lty="dashed",col="red")
mean(sample_cut[[3]])



### SMI model using CN nested MCMC ####

eta=0.75
set.seed(2011320)
sample_smi=CN_smi_sampler(N1=3000,N2=30, log_pr_th=log_prior_th, log_pr_phi=log_prior_phi, l_Z=l_Z, l_Y=l_Y,
                          th0=NULL,phi0=NULL,prop1=MVN_S_prop,prop2=MVN_S_prop,
                          phi_names=c("phi"),th_names=c("theta"),
                          Z=Z,Y=Y,scale1=0.16,scale2=0.1, eta=eta)

phis_smi=sample_smi[[2]]
thetas_smi=sample_smi[[1]]
hist(thetas_smi[,1],breaks=seq(from=-2,to=2,by=0.05),main=paste0("True mean for smi posterior, learning rate:",eta), xlab="phi")
abline(v=0,lty="dashed",col="red")
hist(phis_smi[,1],breaks=seq(from=-2,to=2,by=0.05),main=paste0("Bias for smi posterior, learning rate:",eta), xlab="theta")
abline(v=1,lty="dashed",col="red")
plot(thetas_smi[,1],type="l")
###burn in 50 for standard size training set
plot(phis_smi[,1],type="l")
#burn in 50 by force of thetas
mean(sample_smi[[3]][-(1:100)])
mean(sample_smi[[4]])
mean(sample_smi[[5]][-(1:100)])

### Quality subsection posterior graph (Sec. 5.3.)
etas=c(0,0.25,0.5,0.75,1)
N1=10000
thetas=data.frame(V1=rep(NA,N1-100),X=NA)
phis=data.frame(V1=rep(NA,N1-100),X=NA)

for (i in 1:length(etas)){
  eta=etas[i]
  if(eta<= 0.1){
    sc=1
  }else{
    sc=0.16
  }
  sample_smi=CN_smi_sampler(N1=N1,N2=30, log_pr_th=log_prior_th, log_pr_phi=log_prior_phi, l_Z=l_Z, l_Y=l_Y,
                            th0=NULL,phi0=NULL,prop1=MVN_S_prop,prop2=MVN_S_prop,
                            phi_names=c("phi"),th_names=c("theta"),
                            Z=Z,Y=Y,scale1=sc,scale2=0.1, eta=eta)
  phi=as.vector(sample_smi[[2]])[[1]]
  phis[,i]=phi[-(1:101)]
  theta=as.vector(sample_smi[[1]])[[1]]
  thetas[,i]=theta[-(1:101)]
  print(eta)
}
colnames(phis)=etas
colnames(thetas)=etas

par(mfrow=c(1,1))
cols=rainbow(length(etas))
for (i in 1:length(etas)){
  d_mean=density(thetas[,i])
  if(i==1){
    plot(d_mean,type="l",col=cols[i],ylim=c(0,3),main="True mean posterior", xlab="phi")
    abline(v=0,lty="dashed",col="red")
  }else{
    points(d_mean,type="l",col=cols[i])
  }
}
legend(-2,3,lty=1,legend=c("0 - cut model", "0.25", "0.5", "0.75", "1 - full Bayesian"),col=cols, title="Learning rate")

for (i in 1:length(etas)){
  d_bias=density(phis[,i])
  if(i==1){
    plot(d_bias,type="l",col=cols[i],ylim=c(0,5),main="Bias posterior", xlab="phi")
    abline(v=1,lty="dashed",col="red")
  }else{
    points(d_bias,type="l",col=cols[i])
  }
}
legend(0.5,5,lty=1,legend=c("0 - cut model", "0.25", "0.5", "0.75", "1 - full Bayesian"),col=cols, title="Learning rate")



##### Quantitative analysis (Sec. 5.4.) #####

### Turning a sample into sample of observations of Z, Y (Sec. 5.4.2.) #### 

FB_sample_Z=c()
#note! When adapting to multivariate parameters change lengths into nrows
for (i in 1:nrow(thetas_FB)){
  FB_sample_Z[i]=gen_Z(thetas_FB[i,1])
}

cut_sample_Z=c()
#note! When adapting to multivariate parameters change lengths into nrows
for (i in 1:nrow(thetas_cut)){
  cut_sample_Z[i]=gen_Z(thetas_cut[i,1])
}  

par(mfrow=c(2,2))

hist(FB_sample_Z,breaks=seq(from=-10,to=10,by=0.2))
hist(cut_sample_Z,breaks=seq(from=-10,to=10,by=0.2))

plot(FB_sample_Z,type="l",main="Z sample from full Bayesian posterior", ylab="Z full Bayesian")
plot(cut_sample_Z,type="l", main="Z sample from cut posterior", ylab="Z cut")

acf(FB_sample_Z,main="Z sample autocorrelation from full Bayesian posterior")
acf(cut_sample_Z, main="Z sample autocorrelation from cut posterior")

#we see values are independent so we feel confident to use these
#without any trimming when calculating the energy score


FB_sample_Y=c()
#note! When adapting to multivariate parameters change lengths into nrows
for (i in 1:nrow(thetas_FB)){
  FB_sample_Y[i]=gen_Y(th=thetas_FB[i,1],phi=phis_FB[i,1])
}

cut_sample_Y=c()
#note! When adapting to multivariate parameters change lengths into nrows
for (i in 1:nrow(thetas_cut)){
  cut_sample_Y[i]=gen_Y(th=thetas_cut[i,1],phi=phis_cut[i,1])
}  

hist(FB_sample_Y,breaks=seq(from=-10,to=10,by=0.2))
hist(cut_sample_Y,breaks=seq(from=-10,to=10,by=0.2))

plot(FB_sample_Y,type="l")
plot(cut_sample_Y,type="l",main="Y sample from cut posterior", ylab="Y cut")

acf(FB_sample_Y)
acf(cut_sample_Y,main="Y sample autocorrelation from cut posterior")


#####

### Score variance check (Sec. 5.4.3.) #######


par(mfrow=c(1,1))

size=1
eta=0
set.seed(2011320)
b=1.9
nY=50
nZ=25
nt=20000
Nruns=30
N1=2000
N2=30

Ztest=rnorm(n=nt,mean=0,sd=2)
Ytest=rnorm(n=nt,mean=0,sd=1)

score_rtrain_Z=c()
score_rtrain_Y=c()

if(eta<= 0.1 && size<=4){
  sc=1
}else{
  sc=0.16
}
for(k in 1:Nruns){
  Z=rnorm(n=size*nZ,mean=0,sd=2)
  Y=rnorm(n=size*nY,mean=1,sd=1)
  sample_smi=CN_smi_sampler(N1=N1,N2=N2, log_pr_th=log_prior_th, log_pr_phi=log_prior_phi, l_Z=l_Z, l_Y=l_Y,
                            th0=NULL,phi0=NULL,prop1=MVN_S_prop,prop2=MVN_S_prop,
                            phi_names=c("phi"),th_names=c("theta"),
                            Z=Z,Y=Y,scale1=sc,scale2=0.1, eta=eta)
  
  phis_smi=sample_smi[[2]]
  thetas_smi=sample_smi[[1]]
  phis_smi=phis_smi[-(1:100),]
  thetas_smi=thetas_smi[-(1:100),]
  smi_sample_Y=c()
  #note! When adapting to multivariate parameters change lengths into nrows
  for (z in 1:length(thetas_smi)){
    smi_sample_Y[z]=gen_Y(th=thetas_smi[z],phi=phis_smi[z])
  }
  smi_sample_Z=c()
  #note! When adapting to multivariate parameters change lengths into nrows
  for (z in 1:length(thetas_smi)){
    smi_sample_Z[z]=gen_Z(thetas_smi[z])
  }
  score_rtrain_Z[k]=energy_1D(sample=smi_sample_Z,y=Ztest,b=b)
  score_rtrain_Y[k]=energy_1D(sample=smi_sample_Y,y=Ytest,b=b)
  print(k)
}

par(mfrow=c(2,2))
boxplot(score_rtrain_Z/20000,ylim=c(0,9),main="score of Z, cut")
boxplot(score_rtrain_Y/20000,ylim=c(0,7), main="score of Y, cut")


sqrt(var(score_rtrain_Z/20000))
sqrt(var(score_rtrain_Y/20000))
sqrt(var(score_rtrain_Z/20000)/100)
sqrt(var(score_rtrain_Y/20000)/100)
mean(score_rtrain_Z/20000)
mean(score_rtrain_Y/20000)




######### SMI SCORE CHECK (5.4.4.) ####




#3d grid, consider learning set size x etas

#The table below was computed in parts, thus the code responsible for attaching
#          the scores to already computed matrices for smaller etas.
set.seed(2011320)
b=1.9
sizes=c(0.4, 1, 2, 4, 8)
etas=seq(0.7,0.75,0.05)
oldeta=0.75
nY=50
nZ=25
nt=20000
Ztest=rnorm(n=nt,mean=0,sd=2)
Ytest=rnorm(n=nt,mean=0,sd=1)
Nruns=100
N1=2000
N2=30
score_Z=matrix(0,nrow=length(sizes),ncol=length(etas))
score_Y=matrix(0,nrow=length(sizes),ncol=length(etas))
old_nmZ=paste0("score_table_biased_data/scoreZ_",b,"_to_",oldeta,".csv")
old_nmY=paste0("score_table_biased_data/scoreY_",b,"_to_",oldeta,".csv")
nmZ=paste0("score_table_biased_data/scoreZ_",b,"_to_",tail(etas,1),".csv")
nmY=paste0("score_table_biased_data/scoreY_",b,"_to_",tail(etas,1),".csv")

for (j in 1:length(sizes)){
  size=sizes[j]
  
  
  for(i in 1:length(etas)){
    eta=etas[i]
    if(eta<= 0.1 && size<=4){
      sc=1
    }else{
      sc=0.16
    }
    for(k in 1:Nruns){
      Z=rnorm(n=size*nZ,mean=0,sd=2)
      Y=rnorm(n=size*nY,mean=1,sd=1)
      sample_smi=CN_smi_sampler(N1=N1,N2=N2, log_pr_th=log_prior_th, log_pr_phi=log_prior_phi, l_Z=l_Z, l_Y=l_Y,
                                th0=NULL,phi0=NULL,prop1=MVN_S_prop,prop2=MVN_S_prop,
                                phi_names=c("phi"),th_names=c("theta"),
                                Z=Z,Y=Y,scale1=sc,scale2=0.1, eta=eta)
      
      phis_smi=sample_smi[[2]]
      thetas_smi=sample_smi[[1]]
      phis_smi=phis_smi[-(1:100),]
      thetas_smi=thetas_smi[-(1:100),]
      smi_sample_Y=c()
      #note! When adapting to multivariate parameters change lengths into nrows
      for (z in 1:length(thetas_smi)){
        smi_sample_Y[z]=gen_Y(th=thetas_smi[z],phi=phis_smi[z])
      }
      smi_sample_Z=c()
      #note! When adapting to multivariate parameters change lengths into nrows
      for (z in 1:length(thetas_smi)){
        smi_sample_Z[z]=gen_Z(thetas_smi[z])
      }
      score_Z[j,i]=score_Z[j,i]+energy_1D(sample=smi_sample_Z,y=Ztest,b=b)
      score_Y[j,i]=score_Y[j,i]+energy_1D(sample=smi_sample_Y,y=Ytest,b=b)
      if (k==Nruns/2){
        print("Halfway!") #to monitor progress of the loop
      }
    }
    print(100*((j-1)*length(etas)+i)/(length(etas)*length(sizes))) #to monitor progress of the loop
  }
  #S_Y=score_Y
  #S_Z=score_Z
  #scY=read.csv(old_nmY)[,-1]
  #scZ=read.csv(old_nmZ)[,-1]
  #S_Y=cbind(scY,S_Y)
  #S_Z=cbind(scZ,S_Z)
  #write.csv(S_Y,file=nmY)
  #write.csv(S_Z,file=nmZ)
}

S_Y=score_Y
S_Z=score_Z

scY=scY[,-(1)]
scZ=scZ[,-(1)]
S_Y=cbind(scY,S_Y)
S_Z=cbind(scZ,S_Z)



nmZ=paste0("score_table_biased_data/scoreZ_",b,"_to_",1,".csv")
nmY=paste0("score_table_biased_data/scoreY_",b,"_to_",1,".csv")
colnames(S_Y)=seq(0,1,0.05)
colnames(S_Z)=seq(0,1,0.05)
write.csv(S_Y,file=nmY)
write.csv(S_Z,file=nmZ)


scY=read.csv(nmY)[,-1]
scZ=read.csv(nmZ)[,-1]


S_Y=scY
S_Z=scZ
S_Z=S_Z/(20000*100)
S_Y=S_Y/(20000*100)

etas=seq(0,1,0.05)
sizes=c(0.4, 1, 2, 4, 8)

#Score Y plot
cols=rainbow(5)
for (i in 1:length(sizes)){
  if(i==1){
    plot(etas,S_Y[i,],type="l",col=cols[i],ylim=c(2,3.5),ylab="Score",main="Score of Y")
  }else{
    points(etas,S_Y[i,],type="l",col=cols[i])
  }
}
legend(0.75,2.5,lty=1,legend=c("10 Z, 20 Y", "25 Z, 50 Y","50 Z, 100Y", "100Z, 200Y", "200 Z, 400 Y"),col=cols, title="Train set size")

#Score Z plot
cols=rainbow(5)
for (i in 1:length(sizes)){
  if(i==1){
    plot(etas,S_Z[i,],type="l",col=cols[i],ylim=c(7,8.4),ylab="Score",main="Score of Z")
  }else{
    points(etas,S_Z[i,],type="l",col=cols[i])
  }
}
legend(0.75,7.5,lty=1,legend=c("10 Z, 20 Y", "25 Z, 50 Y","50 Z, 100Y", "100Z, 200Y", "200 Z, 400 Y"),col=cols, title="Train set size")


