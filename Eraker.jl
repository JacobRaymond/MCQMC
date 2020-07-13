using CSV, Distributions, LinearAlgebra, Plots

#Import the data  (https://fred.stlouisfed.org/graph/?g=r6hR)
Y_dat=CSV.read("/Users/JacobRaymond 1/Library/Mobile Documents/com~apple~CloudDocs/Maitrise/Papier/MCQMC/WTB6MS.csv").WTB6MS
#Y=Y./100
Y_dat=Y_dat./100
Y_dat=Y_dat[(521-(2*52)):520]

#Set variables
m=5
del=1/m

#Simulate initial parameter values
sig=rand()
sig_vec=[sig]
a=rand(Uniform(0, sig))
b=rand(Uniform(0, sig))
phi=[[a,b]]

#Initial path
Y=[]
for i in 1:(length(Y_dat)-1)
    Y_iter=[Y_dat[i]]
    for k in 1:(m-1)
        push!(Y_iter, Y_dat[i]+((Y_dat[i+1]-Y_dat[i])/m)*k)
    end
    push!(Y, Y_iter)
end
Y=collect(Iterators.flatten(Y))
push!(Y, last(Y))

#Path length
T=length(Y)

#### MCMC ####

#Number of iterations
N=2039

#Gibbs Sampler
for k in 1:N
    global Y
    global a
    global b
    global del
    global sig

    #Augmented data
    Y_new=[Y[1]]
    for i in 2:(T-1)
        mean_y=(Y[i-1]+Y[i+1])/2
        sig_y=sqrt(0.5*sig*del)

        prop_y=rand(Normal(Y[i], var(Y)))

        yratio=pdf(Normal(mean_y, sig_y), prop_y)/pdf(Normal(mean_y, sig_y), Y[i])
        yratio=pdf(Normal(Y[i], var(Y)),  Y[i])/pdf(Normal(Y[i], var(Y)), prop_y)*yratio

        if yratio>rand()
            push!(Y_new, prop_y)
        else
            push!(Y_new, Y[i])
        end
    end
    push!(Y_new, last(Y))
    Y=Y_new

    #Calculate some variables for the distribution of (a,b)
    y=diff(Y)./sqrt.(Y[1:(T-1)]*del)
    X=hcat((Y[1:length(Y)-1]).^(-0.5), (Y[1:length(Y)-1]).^(0.5))
    X=sqrt(del).*X

    #Mean and variance
    var_ab=inv(Transpose(X)*X)
    mu=var_ab*Transpose(X)*y

    #Generate new value
    phi_it=rand(MvNormal(mu, sig*var_ab))


    #Save values
    a=phi_it[1]
    b=phi_it[2]
    push!(phi, [a,b])

    #Calculate some variables for the distribution of sigma^2
    E=T-2
    F=mean((y.-(X*mu)).^2)

    #Generate a new value of sigma^2
    sig=rand(InverseGamma(E,F))
    push!(sig_vec, sig)
end

#Results

println("MCMC Results:")

#Estimate for a and b (500 observation burn)
println("(a, b)=", mean(phi[500:N]))

#Estimate for sigma^2
println("sigma^2=", mean(sqrt.(sig_vec[500:N])))

#Estimate for the mean of the Y
a_es=mean(phi[500:N])[1]
b_es=mean(phi[500:N])[2]

#Ratios (should be close)
println("Ratio:", a_es/b_es)
println(mean(Y))

#### MCQMC ####
include("LCG.jl")
include("MVN_QMC.jl")

#Generate points
N=2039
mult=995
u=lcg(N, mult, 3)

#Gibbs Sampler
for k in 1:N
    global Y
    global a
    global b
    global del
    global sig

    #Augmented data
    Y_new=[Y[1]]
    for i in 2:(T-1)
        mean_y=(Y[i-1]+Y[i+1])/2
        sig_y=sqrt(0.5*sig*del)

        prop_y=rand(Normal(Y[i], var(Y)))

        yratio=pdf(Normal(mean_y, sig_y), prop_y)/pdf(Normal(mean_y, sig_y), Y[i])
        yratio=pdf(Normal(Y[i], var(Y)),  Y[i])/pdf(Normal(Y[i], var(Y)), prop_y)*yratio

        if yratio>rand()
            push!(Y_new, prop_y)
        else
            push!(Y_new, Y[i])
        end
    end
    push!(Y_new, last(Y))
    Y=Y_new

    #Calculate some variables for the distribution of (a,b)
    y=diff(Y)./sqrt.(Y[1:(T-1)]*del)
    X=hcat((Y[1:length(Y)-1]).^(-0.5), (Y[1:length(Y)-1]).^(0.5))
    X=sqrt(del).*X

    #Mean and variance
    var_ab=inv(Transpose(X)*X)
    mu=var_ab*Transpose(X)*y

    #Generate new value
    phi_it=rand(MvNormal(mu, sig*var_ab))
    phi_it=MVN_QMC(u[k][1:2], mu, sig*var_ab)

    #Save values
    a=phi_it[1]
    b=phi_it[2]
    push!(phi, [a,b])

    #Calculate some variables for the distribution of sigma^2
    E=T-2
    F=mean((y.-(X*mu)).^2)

    #Generate a new value of sigma^2
    sig=quantile(InverseGamma(E,F), u[k][3])
    push!(sig_vec, sig)
end

#Results

println("MC-QMC Results:")

#Estimate for a and b (500 observation burn)
println("(a, b)=", mean(phi[500:1500]))

#Estimate for sigma^2
println("sigma^2=", mean(sqrt.(sig_vec[500:1500])))

#Estimate for the mean of the Y
a_es=mean(phi[500:1500])[1]
b_es=mean(phi[500:1500])[2]

#Ratios (should be close)
println("Ratio:", a_es/b_es)
println(mean(Y))
