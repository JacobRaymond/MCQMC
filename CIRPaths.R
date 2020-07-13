library(tidyverse)

#Data
Y=read.csv("WTB6MS.csv")%>%
  slice(tail(row_number(), 104)) %>% #We select the last two years
  select(WTB6MS) %>%
  mutate(WTB6MS=WTB6MS/100)%>%
  as_vector()

#Model Parameters
a=0.00094 
b=-0.474010
sigma=0.0002253

#Confidence
alpha=0.05


#Simulate path
pathsim=function(a, b, sigma_sq,  Y, delta=1, plength=length(Y)){
  
  #Vector to house observations
  p=vector(length=plength)
  
  #Initial value
  p[1]=Y[1]
  
  #Simulate path using CIR Discretization
  for(i in 2:plength){
    p[i]<--1
    while(p[i]<0){
      p[i]<-p[i-1]+(a+b*p[i-1])*delta+sqrt(sigma_sq*delta*p[i-1])*rnorm(1)
    }
  }
  p
}


CIRPaths=function(a, b, sigma_sq, Y, paths=100, delta=1,  plength=length(Y), alpha=0.05){
  
  #Simulate multiple paths and extract the confidence interval bounds
  pathlist= map(seq_len(paths), ~pathsim(a, b, sigma, Y, delta, plength)) %>% 
    bind_cols() %>% 
    mutate(lb=apply(., 1, quantile, probs = alpha), ub=apply(., 1, quantile, probs = 1-alpha))
  
  #Plot output
  ggplot(data=pathlist, aes(x=seq_len(plength)))+
    geom_line(data=NULL, aes(x=seq_len(length(Y)), y=Y, colour="Y"))+
    geom_ribbon(aes( ymin = lb, ymax = ub, fill="Confidence Band"),  alpha = .25)+
    labs(x="", y="", fill="")+
    scale_y_continuous(labels = scales::percent)+
    scale_colour_manual("", values="black")+
    theme(legend.position = "bottom")
}

CIRPaths(a,b,sigma, Y, alpha=0.1)

