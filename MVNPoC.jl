using Distributions, PrettyTables, LinearAlgebra
include("LCG.jl")
include("MVN_QMC.jl")

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

    #MSE (the true value is 0)
    MSE=mean(map(x->x.^2, sim_res))

    (Mean, MSE)
end


#### QMC MCMC ####

function qmc_mcmc(prop)
    sim_res=[]

    #Correlation Matrix
    Sigma=Matrix(Diagonal(repeat([0.1], 20)))

    for j in 1:100

        #Initial point - simulate a random point between [0,1]
        x=rand(20)
        mh_pts=[]

        #Generate rqmc points
        u=lcg(1021, 65, 21)

        for i in 1:1021

            #Simulate new point
            y=MVN_QMC(u[i][1:20], x, Sigma)

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

    #MSE (the true value is 0)
    MSE=mean(map(x->x.^2, sim_res))

    (Mean, MSE)
end

#### Simulation Study ###

#Random Walk
function prop_rw(w)
    MvNormal(w, Diagonal(repeat([0.1], 20)))
end
res_mc=crude_mcmc(prop_rw)
res_qmc=qmc_mcmc(prop_rw)

#Output in table
data=Any[round.(res_mc[1], digits=3) round.(res_qmc[1], digits=3) round.(res_mc[2], digits=3) round.(res_qmc[2], digits=3) round.(res_mc[2]./res_qmc[2], digits=3)]
header=["Mean (MC)" "Mean(QMC)" "MSE (MC)" "MSE (QMC)" "Ratio"]
pretty_table(data, header)
