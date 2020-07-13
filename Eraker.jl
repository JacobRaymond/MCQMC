using CSV, Distributions, LinearAlgebra, PrettyTables

#Import the data  (https://fred.stlouisfed.org/graph/?g=r6hR)
Y_dat=CSV.read("/Users/JacobRaymond 1/Library/Mobile Documents/com~apple~CloudDocs/Maitrise/Papier/MCQMC/WTB6MS.csv").WTB6MS
#Y=Y./100
Y_dat=Y_dat./100
Y_dat=Y_dat[(521-(2*52)):520]

#Set variables
m=5
del=1/m

#### MCMC ####

#Vector to save outputs
phi_gibbs_mc=[]
sig_gibbs_mc=[]

for l in 1:100
    global del

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

    #Number of iterations
    N=1021

    #Gibbs Sampler
    for k in 1:N

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

    #Save
    push!(phi_gibbs_mc, mean(phi))
    push!(sig_gibbs_mc, mean(sig_vec))
end


#### MC-QMC ####
include("LCG.jl")
include("MVN_QMC.jl")

#Vector to save outputs
phi_gibbs_qmc=[]
sig_gibbs_qmc=[]

for l in 1:100
    global del

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

    #Generate points
    N=1021
    mult=65
    u=lcg(N, mult, 3)

    #Gibbs Sampler
    for k in 1:N

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

    #Save
    push!(phi_gibbs_qmc, mean(phi))
    push!(sig_gibbs_qmc, mean(sig_vec))
end


#### Results ####

#Show Results
data_res_mc=vcat(mean(map(x->x[1], phi_gibbs_mc)), mean(map(x->x[2], phi_gibbs_mc)), mean(sig_gibbs_mc))
data_res_qmc=vcat(mean(map(x->x[1], phi_gibbs_qmc)), mean(map(x->x[2], phi_gibbs_qmc)), mean(sig_gibbs_qmc))
cols=["a", "b", "σ^2"]
data_res=hcat(cols, round.(data_res_mc, digits=7), round.(data_res_qmc, digits=7))
header_res=["Parameter","MCMC", "MC-QMC"]
pretty_table(data_res, header_res)

#Variance Reduction
data_var_mc=vcat(var(map(x->x[1], phi_gibbs_mc)), var(map(x->x[2], phi_gibbs_mc)), var(sig_gibbs_mc))
data_var_qmc=vcat(var(map(x->x[1], phi_gibbs_qmc)), var(map(x->x[2], phi_gibbs_qmc)), var(sig_gibbs_qmc))
cols=["a", "b", "σ^2"]
ratio=data_var_mc./data_var_qmc
data_res=hcat(cols, round.(data_var_mc, digits=10), round.(data_var_qmc, digits=10), round.(ratio, digits=3))
header_res=["Parameter","Variance (MCMC)", "Variance (MC-QMC)", "Ratio"]
pretty_table(data_res, header_res)
