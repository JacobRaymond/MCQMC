using Distributions, PrettyTables, LinearAlgebra
include("LCG.jl")
include("IW_QMC.jl")

#### Crude MCMC ####

        #Parameter for simulation study; unnformative prior
        v=5
        X=Matrix(v*Diagonal(ones(3)))

        #Simulate data
        df=rand(MvNormal(Diagonal(ones(3))), 100)

        #Initial point - simulate a random matrix and make it hermitian
        x=rand(3,3)
        x=x*transpose(x)
        mh_pts=[x]

        for i in 1:32749
            global x

            #Simulate new point
            y=rand(InverseWishart(v, X))

            #Calculate Metropolis ratio
            mh=((sum(logpdf(MvNormal(y), df)))*pdf(InverseWishart(v, X), y))/((sum(logpdf(MvNormal(x), df)))*pdf(InverseWishart(v, X), x))

            #Calculate acceptance probability
            x=ifelse(mh>rand(), y, x)

            #Save value
            push!(mh_pts, x)
        end

        println(mean(mh_pts))
    
