using Distributions, PrettyTables, LinearAlgebra
include("LCG.jl")
include("IW_QMC.jl")

#### Crude MCMC ####
sim_res=[]

for k in 1:100

        #Parameter for simulation study; unnformative prior
        v=3
        X=Matrix(v*Diagonal(ones(3)))

        #Simulate data
        df=rand(MvNormal(Diagonal(ones(3))), 100)

        #Initial point - simulate a random matrix and make it hermitian
        x=rand(3,3)
        x=x*transpose(x)
        mh_pts=[x]

        for i in 1:32749

            #Simulate new point
            y=rand(InverseWishart(v, X))

            #Calculate Metropolis ratio
            mh=((sum(logpdf(MvNormal(y), df)))*pdf(InverseWishart(v, X), y))/((sum(logpdf(MvNormal(x), df)))*pdf(InverseWishart(v, X), x))

            #Calculate acceptance probability
            x=ifelse(mh>rand(), y, x)

            #Save value
            push!(mh_pts, x)
        end

    push!(sim_res, mean(mh_pts))
end

#Mean (Should be uniform)
println("Mean Estimate (MC)")
pretty_table(round.(mean(sim_res), digits=5), tf = borderless, noheader = true)

#MSE
MSE_mat=map(x->x-Diagonal(ones(3)), sim_res)
MSE_mat_mc=map(x->x.^2, MSE_mat)
println("MSE (MC)")
pretty_table(round.(mean(MSE_mat_mc), digits=5), tf = borderless, noheader = true)


#### MCQMC ####
sim_res=[]

for k in 1:100

        #Parameter for simulation study; unnformative prior
        v=3
        X=Matrix(v*Diagonal(ones(3)))

        #Simulate data
        df=rand(MvNormal(Diagonal(ones(3))), 100)

        #Initial point - simulate a random matrix and make it hermitian
        x=rand(3, 3)
        x=x*transpose(x)
        mh_pts=[x]

        #QMC Points
        u=lcg(32749, 219, 7)

        for i in 1:32749

            #Simulate new point
            y=IW_QMC(u[i][1:6], X, v)


            #Calculate Metropolis ratio
            mh=((sum(logpdf(MvNormal(y), df)))*pdf(InverseWishart(v, X), y))/((sum(logpdf(MvNormal(x), df)))*pdf(InverseWishart(v, X), x))

            #Calculate acceptance probability
            x=ifelse(mh>u[i][7], y, x)

            #Save value
            push!(mh_pts, x)
        end

    push!(sim_res, mean(mh_pts))
end

#Mean (Should be uniform)
println("Mean Estimate (QMC)")
pretty_table(round.(mean(sim_res), digits=5), tf = borderless, noheader = true)

#MSE
MSE_mat=map(x->x-Diagonal(ones(3)), sim_res)
MSE_mat_qmc=map(x->x.^2, MSE_mat)
println("MSE (QMC)")
pretty_table(round.(mean(MSE_mat_qmc), digits=5), tf = borderless, noheader = true)

####Variance Reduction####
println("Variance Reduction:")
MSE_red=MSE_mat_mc ./ MSE_mat_qmc
pretty_table(round.(abs.(mean(MSE_red)), digits=5), tf = borderless, noheader = true)
