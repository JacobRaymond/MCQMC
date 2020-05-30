using CSV, Distributions, LinearAlgebra

#Import the data  (https://fred.stlouisfed.org/graph/?g=r6hR)
Y_dat=CSV.read("/Users/JacobRaymond 1/Desktop/WTB6MS.csv")
Y_dat=Y_dat.WTB6MS

#Set variables
M=20
del_plus=1/52
del=del_plus/(M+1)
T=length(Y_dat)-1

#Simulate initial parameter values
a=rand()
b=rand()
phi=[[a,b]]
sig=var(Y_dat) #Initial value: the variance of the yields
sig_vec=[sig]

#Create Y
Y=Array[];
for i in 1:length(Y_dat)
    push!(Y, [Y_dat[i]])
end
pop!(Y)

#Function to simulate augmented data
function augsim(y, a, b, del, sig)
    mean_sim=y+(a-b*y)*del
    var_sim=sig*del*y
    rand(Normal(mean_sim, sqrt(var_sim)))
end

#Number of iterations
N=1500

#Gibbs Sampler
for k in 1:N
    global augsim
    global Y
    global a
    global b
    global del
    global sig

    #Initialize
    Y_star=deepcopy(Y);
    #Simulate augmented data
    for i in 1:T
        for j in 1:M
            push!(Y_star[i], augsim(Y_star[i][j], a, b, del, sig))
        end
    end

    #Calculate some variables for the distribution of (a,b)
    A=sum(map(y->sum( map(x->1/x, y)), Y_star))
    B=sum(map(x->sum(x), Y_star))
    C=sum(map(x->sum(diff(x)./x[1:M]), Y_star))
    D=sum(map(x->sum(-diff(x)), Y_star))

    #Create matrix components
    a11=(del*A)/sig
    a12=-(del*(T-1)*(M+1))/sig
    a22=(del*B)/sig

    #Find the mean
    mu=[(a22*C-a12*D)/(a11*a22-a12^2), (a11*D-a12*C)/(a11*a22-a12^2)]

    #Variance-covariance matrix
    vcov= (1/(a11*a22-a12^2))*[a22 -a12; -a12 a11]
    phi_it=rand(MvNormal(mu, vcov))

    #Save values
    a=phi_it[1]
    b=phi_it[2]
    push!(phi, phi_it)

    #Calculate some variables for the distribution of sigma^2
    E=0.5*(T-1)*(M+1)
    function F_fun(x)
        val_tot=0
        for i in 1:M
            val_tot+=((x[i+1]-(a+b*x[i])*del)^2)/(2*x[i])
        end
        val_tot
    end
    F=sum(map(x->F_fun(x), Y_star))

    #Generate a new value of sigma^2
    sig=rand(InverseGamma(E, F))
    push!(sig_vec, sig)
end

#Estimate for a and b (500 observation burn)
mean(phi[500:1500])

#Estimate for sigma^2
mean(sig_vec[500:1500])
