using Distributions, PrettyTables, LinearAlgebra
include("LCG.jl")

#### Crude MCMC ####

function crude_mcmc(prop)
    sim_res=[]
    for j in 1:100

        #Initial point - simulate a random point between [0,1]
        x=rand(20)
        mh_pts=[]

        for i in 1:1021

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

    #Variance
    Variance=var(sim_res)

    (Mean, Variance)
end


#### QMC MCMC ####

function qmc_mcmc(prop)
    sim_res=[]

    #Cholesky decomposition
    Sigma=Matrix(Diagonal(repeat([0.1], 20)))
    A=cholesky(Sigma)
    A=A.U

    for j in 1:100

        #Initial point - simulate a random point between [0,1]
        x=rand(20)
        mh_pts=[]

        #Generate rqmc points
        u=lcg(1021, 65, 21)

        for i in 1:1021

            #Simulate new point
            z=map(x->quantile(Normal(), x), u[i][1:20])
            y=x+A*z

            #Calculate Metropolis ratio
            mh=(pdf(MvNormal(Diagonal(ones(20))), y)*pdf(prop(y), x))/(pdf(MvNormal(Diagonal(ones(20))), x)*pdf(prop(x), y))

            #Calculate acceptance probability
            x=ifelse(mh>u[i][21], y, x)

            #Save value
            push!(mh_pts, x)
        end

        #Calculate Mean Estimator
        push!(sim_res, mean(mh_pts))
    end

    #Mean
    Mean= mean(sim_res)

    #MSE (true value is 0)
    Variance=var(sim_res)

    (Mean, Variance)
end

function prop_rw(w)
    MvNormal(w, Diagonal(repeat([0.1], 20)))
end

res_mc=crude_mcmc(prop_rw)
res_qmc=qmc_mcmc(prop_rw)
println(res_mc.-res_qmc)
