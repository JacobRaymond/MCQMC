using Distributions, PrettyTables
include("LCG.jl")

#### Crude MCMC ####

function crude_mcmc(prop)
    sim_res=Float64[]
    for j in 1:300

        #Initial point - simulate a random point between [0,1]
        x=rand()
        mh_pts=Float64[]

        for i in 1:65521

            #Simulate new point
            y=rand(prop(x))

            #Calculate Metropolis ratio
            mh=(pdf(Normal(), y)*pdf(prop(y), x))/(pdf(Normal(), x)*pdf(prop(x), y))

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

    #MSE (true value is 0)
    MSE= mean(sim_res.^2)

    (Mean, MSE)
end


#### QMC MCMC ####

function qmc_mcmc(prop)
    sim_res=Float64[]
    for j in 1:300

        #Create a CUD sequence
        rqmc=lcg(65521, 17364, 2)

        #Extract proposals
        y=map(w->w[1], rqmc)

        #Extract acceptance probabilities
        u=map(w->w[2], rqmc)


        #Initial point - simulate a random point between [0,1]
        x=rand()
        mh_pts=Float64[]

        for i in 1:65521

            #Simulate new point
            y_prop=quantile(prop(x), y[i])

            #Calculate Metropolis ratio
            mh=(pdf(Normal(), y_prop)*pdf(prop(y_prop), x))/(pdf(Normal(), x)*pdf(prop(x), y_prop))

            #Calculate acceptance probability
            x=ifelse(mh>u[i], y_prop, x)

            #Save value
            push!(mh_pts, x)
        end

        #Calculate Mean Estimator
        push!(sim_res, mean(mh_pts))
    end

    #Mean
    Mean= mean(sim_res)

    #MSE (true value is 0)
    MSE= mean(sim_res.^2)

    (Mean, MSE)
end



#### Simulation Study ####

#Independent proposal
function prop_ind(w)
    Normal(0, 2.4)
end
res_ind_crude=crude_mcmc(prop_ind)
res_ind_qmc=qmc_mcmc(prop_ind)

#Random walk
function prop_rw(w)
    Normal(w, 2.4)
end
res_rw_crude=crude_mcmc(prop_rw)
res_rw_qmc=qmc_mcmc(prop_rw)


#Output results
data=Any["Independent" "Crude MC" round.(res_ind_crude, digits=6)[1] round.(res_ind_crude, digits=6)[2] ;
"Independent" "QMC" round.(res_ind_qmc, digits=6)[1] round.(res_ind_qmc, digits=6)[2] ;
"Random Walk" "Crude MC" round.(res_rw_crude, digits=6)[1] round.(res_rw_crude, digits=6)[2] ;
"Random Walk" "QMC" round.(res_rw_qmc, digits=6)[1] round.(res_rw_qmc, digits=6)[2] ;
]
header=["Proposal" "MCMC Type" "Mean" "MSE"]
pretty_table(data, header)

#Improvement in MSE
ind_mse=round(res_ind_crude[2]/res_ind_qmc[2], digits=2)
rw_mse=round(res_rw_crude[2]/res_rw_qmc[2], digits=2)
println("Using quasi-Monte Carlo inputs reduced the MSE by a factor of ", ind_mse, " for the independence sampler, and by ",
rw_mse, " for the random walk example.")
