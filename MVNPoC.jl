using Distributions, PrettyTables, LinearAlgebra
include("LCG.jl")

#### Crude MCMC ####

function crude_mcmc(prop)
    sim_res=[]
    for j in 1:100

        #Initial point - simulate a random point between [0,1]
        x=rand(20)
        mh_pts=[]

        for i in 1:65521

            #Simulate new point
            y=rand(prop(x))

            #Calculate Metropolis ratio
            mh=(pdf(MvNormal(Diagonal(ones(20))), y)*pdf(prop(y), x))/(pdf(MvNormal(Diagonal(ones(20))), x)*pdf(prop(x), y))

            #Calculate acceptance probability
            x=ifelse(mh>rand(), y, x)

            #Save value
            push!(mh_pts, x)
        end

        #Calculate Mean Estimator
        push!(sim_res, mean(mh_pts))
    end

    #Mean
    Mean= mean(sim_res)
    #Var=mean(var(sim_res))
    #(Mean, Var)

    #MSE (true value is 0)
    MSE= map(x->mean(x.^2), sim_res)

    (Mean, MSE)
end

function prop_rw(w)
    MvNormal(w, Diagonal(repeat([0.1], 20)))
end

res_mc=crude_mcmc(prop_rw)
