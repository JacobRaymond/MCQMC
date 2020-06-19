using CSV, Distributions, LinearAlgebra

#Import the data  (https://fred.stlouisfed.org/graph/?g=r6hR)
Y_dat=CSV.read("/Users/JacobRaymond 1/Desktop/WTB6MS.csv").WTB6MS./100


#Set variables
M=20
del_plus=1/52
del=del_plus/(M+1)
T=length(Y_dat)

#Simulate initial parameter values
sig=0.0014 #Initial value: the variance of the yields
sig_vec=[sig]
a=rand(Uniform(0, sqrt(sig)))
b=rand(Uniform(0, sqrt(sig)))
phi=[[a,b]]

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
    candidate=rand(Normal(mean_sim, sqrt(var_sim)))
    if candidate>0
        candidate
    else y end
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
    if M>0
        for i in 1:T-1
            for j in 1:M
                push!(Y_star[i], augsim(Y_star[i][j], a, b, del, sig))
            end
        end
    end
    Y_star=collect(Iterators.flatten(Y_star))

    #Calculate some variables for the distribution of (a,b)
    A=sum(1 ./Y_star)
    B=sum(Y_star)
    C=sum(diff(Y_star)./Y_star[1:length(Y_star)-1])
    D=sum(-diff(Y_star))

    #Create matrix components
    a11=(del*A)/sig
    a12= -(del*(T-1)*(M+1))/sig
    a22=(del*B)/sig

    #Find the mean
    mu=[(a22*C-a12*D)/(a11*a22-a12^2), (a11*D-a12*C)/(a11*a22-a12^2)]

    #Variance-covariance matrix
    vcov= [a11 a12; a12 a22]
    vcov=inv(vcov)
    phi_it=rand(MvNormal(mu, vcov))


    #Save values
    a=phi_it[1]
    b=phi_it[2]
    push!(phi, [a,b])

    #Calculate some variables for the distribution of sigma^2
    E=0.5*(T-1)*(M+1)

    F=sum(((diff(Y_star).-(a.-b*Y_star[1:length(Y_star)-1])*del).^2)./(2*Y_star[1:length(Y_star)-1]))

    #Generate a new value of sigma^2
    sig=rand(InverseGamma(E, F))
    push!(sig_vec, sig)
end

#Estimate for a and b (500 observation burn)
println(mean(phi[500:1500]))

#Estimate for sigma^2
println(mean(sig_vec[500:1500]))

#Estimate for the mean of the y_dat
a_es=mean(phi[500:1500])[1]
b_es=mean(phi[500:1500])[2]
print(a_es/b_es) #Should be close to 0.023
