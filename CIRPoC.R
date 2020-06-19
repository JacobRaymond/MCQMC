library(mvtnorm)
library(invgamma)

a_it=c()
b_it=c()
sig_it=c()

#Import data (last two years of the data)
y_dat=read.csv("/Users/JacobRaymond 1/Library/Mobile Documents/com~apple~CloudDocs/Maitrise/Papier/MCQMC/WTB6MS.csv")
y_dat=y_dat$WTB6MS/100
y_dat=y_dat[417:520]

#Set variables
M=20
del_plus=1/52
del=del_plus/(M+1)
T=length(y_dat)-1

#### Method 1: Acceptance Rejection ####
start<-Sys.time()
for(l in 1:10){
  #Starting values
  sig=sd(y_dat)
  sig_gibbs<-sig
  a=runif(max = sig, 1)
  b=runif(max=sig, 1)
  phi_gibbs=c(a,b)
  
  for(k in 1:1500){
    #Generate intermediary values
    y_star=as.list(y_dat)
    for(j in 1:T){
      y<-y_dat[[j]]
      for(i in 1:M){
        new_y=-1
        while(new_y<0){
          new_y=rnorm(mean = y[i]+(a-b*y[i])*del, sd = sig*sqrt(del*y[i]), n=1)
        }
        y<-c(y, new_y)
      }
      y_star[[j]]<-y
    }
    y_star<-unlist(y_star)
    
    #Calculate values for distribution of a and b
    A=sum(1/y_star)
    B=sum(y_star)
    C=-sum(1-y_star[-1]/y_star[-length(y_star)])
    D=sum(-diff(y_star))
    a11=del*A/sig^2
    a12=-(del/sig^2)*(T-1)*(M+1)
    a22=del*B/sig^2
    mu1=(a22*C-a12*D)/(a11*a22-a12^2)
    mu2=(-a12*C+a11*D)/(a11*a22-a12^2)
    
    #Mean
    mu=c(mu1, mu2)
    
    #Variance
    Lambda=matrix(c(a11, a12, a12, a22), nrow = 2, ncol=2)
    Lambda=solve(Lambda)
    
    #Generate new value for a, b
    phi=rmvnorm(n=1, mean = mu, sigma = Lambda)
    phi_gibbs<-rbind(phi_gibbs, phi)
    a=phi[1]
    b=phi[2]
    
    #Calculate values for distribution of sig
    E=(T-1)*(M+1)/2
    numerator=diff(y_star)-(a-b*y_star[-length(y_star)])*del
    F=sum((numerator^2)/(2*y_star[-length(y_star)]))
    sig<-rinvgamma(n=1, E, F)
    sig_gibbs<-c(sig_gibbs,sig)
  }
  
  #Estimate sigma
  sig_it=c(sig_it, mean(sig_gibbs[500:1500]))
  
  #Estimate a
  a_it=c(a_it, mean(phi_gibbs[500:1500,1]))
  
  #Estimate b
  b_it=c(b_it, mean(phi_gibbs[500:1500,2]))
}
end<-Sys.time()
end-start

#Estimate of a (1.583 X 10^-3)
mean(a_it)

#Estimate of b (6.944 X 10^-3)
mean(b_it)

#Ratio (0.2618)
mean(a_it)/mean(b_it)

#Estimate of sigma (√0.0398)
mean(sig_it)

#### Method 2: Transformation ####
start<-Sys.time()
for(l in 1:10){
  #Starting values
  sig=sd(y_dat)
  sig_gibbs<-sig
  a=runif(max = sig, 1)
  b=runif(max=sig, 1)
  phi_gibbs=c(a,b)
  
  for(k in 1:1500){
    #Generate intermediary values
    y_star=as.list(y_dat)
    for(j in 1:T){
      y<-y_dat[[j]]
      for(i in 1:M){
          mean_y=y[i]+(a-b*y[i])*del
          sd_y=sig*sqrt(del*y[i])
          
          #Transformation
          transfo=pnorm(0, mean = mean_y, sd=sd_y)+runif(1)*(1-pnorm(0, mean = mean_y, sd=sd_y))
          new_y=qnorm(transfo, mean = mean_y, sd=sd_y)
          
        y<-c(y, new_y)
        }
      y_star[[j]]<-y
    }
    y_star<-unlist(y_star)
    
    #Calculate values for distribution of a and b
    A=sum(1/y_star)
    B=sum(y_star)
    C=-sum(1-y_star[-1]/y_star[-length(y_star)])
    D=sum(-diff(y_star))
    a11=del*A/sig^2
    a12=-(del/sig^2)*(T-1)*(M+1)
    a22=del*B/sig^2
    mu1=(a22*C-a12*D)/(a11*a22-a12^2)
    mu2=(-a12*C+a11*D)/(a11*a22-a12^2)
    
    #Mean
    mu=c(mu1, mu2)
    
    #Variance
    Lambda=matrix(c(a11, a12, a12, a22), nrow = 2, ncol=2)
    Lambda=solve(Lambda)
    
    #Generate new value for a, b
    phi=rmvnorm(n=1, mean = mu, sigma = Lambda)
    phi_gibbs<-rbind(phi_gibbs, phi)
    a=phi[1]
    b=phi[2]
    
    #Calculate values for distribution of sig
    E=(T-1)*(M+1)/2
    numerator=diff(y_star)-(a-b*y_star[-length(y_star)])*del
    F=sum((numerator^2)/(2*y_star[-length(y_star)]))
    sig<-rinvgamma(n=1, E, F)
    sig_gibbs<-c(sig_gibbs,sig)
  }
  
  #Estimate sigma
  sig_it=c(sig_it, mean(sig_gibbs[500:1500]))
  
  #Estimate a
  a_it=c(a_it, mean(phi_gibbs[500:1500,1]))
  
  #Estimate b
  b_it=c(b_it, mean(phi_gibbs[500:1500,2]))
}
end<-Sys.time()
end-start

#Estimate of a (1.583 X 10^-3)
mean(a_it)

#Estimate of b (6.944 X 10^-3)
mean(b_it)

#Ratio (0.2618)
mean(a_it)/mean(b_it)

#Estimate of sigma (√0.0398)
mean(sig_it)
