using Distributions, PrettyTables
include("Korobov.jl")

# Load data arrays
s=[5,1,5,14,3,19,1,1,4,22]
t=[94.320, 15.720, 62.880, 125.760, 5.240, 31.440, 1.048, 1.048, 2.096, 10.480]

#Variables
a=1.802
del=1
gam=0.1


#### Crude Method (see DOI: 10.2307/2289776) ####

#Initial values
init=s./t #MLE
b=(sum(init)+del)/(gam+10*a-1)#Mean of a IG(gamma+10a, ∑lambda+∂)

#Array to house the observations
lambda_gibbs= Array{Float64}(undef, 0, 10)
b_gibbs=Float64[]

#Gibbs Sampler

m=300 #Number of times the sampler is repeated
n=1021 #Number of iterations of the sampler

for k in 1:m
    global lambda_gibbs
    for j in 1:n
        global b
        global lambda
        global s
        global t
        local lambda_gibbs

        #Draw the lambdas
        lambda=Float16[rand(Gamma(a+s[1], 1/(t[1]+1/b)))]
        for i in 2:10
            push!(lambda, rand(Gamma(a+s[i], 1/(t[i]+1/b))))
        end

        #Draw the beta
        b=rand(InverseGamma(gam+10*a, sum(lambda)+del))
    end

    #Save the values
    lambda_gibbs=[lambda_gibbs; lambda']
    push!(b_gibbs, b)
end

#Save outputs
lambda_mc=lambda_gibbs
b_mc=b_gibbs

#### RQMC Method (see  DOI: 10.1073/pnas.0409596102)####

#Initial values
init=s./t #MLE
b=(sum(init)+del)/(gam+10*a-1)#Mean of a IG(gamma+10a, ∑lambda+∂)

#Array to house the observations
lambda_gibbs= Array{Float64}(undef, 0, 10)
b_gibbs=Float16[]

#Draw the lambdas
m=300 #The algorithm is repeated 300 times, as per Owen and Tribble
n=1021 #Number of iterations of the sampler

for k in 1:m
    global lambda_gibbs
    #Generate RMQC points
    rqmc=Korobov(1021, 65, 11)
    for j in 1:n
        global b
        global lambda
        global s
        global t

        #Draw the lambdas
        lambda=Float16[quantile(Gamma(a+s[1], 1/(t[1]+1/b)), rqmc[j][1])]
        for i in 2:10
            push!(lambda, quantile(Gamma(a+s[i], 1/(t[i]+1/b)), rqmc[j][i]))
        end

        #Draw the beta
        b=quantile(InverseGamma(gam+10*a, sum(lambda)+del), rqmc[j][11])
    end

    #Save the values
    lambda_gibbs=[lambda_gibbs; lambda']
    push!(b_gibbs, b)
end

#Save outputs
lambda_qmc=lambda_gibbs
b_qmc=b_gibbs

#### Output ####

#Calculate estimators, MC Method
result_mc=mean.([lambda_mc[:,j] for j in 1:size(lambda_mc,2)])
push!(result_mc, mean(b_mc))
println("The estimators, calculated using the crude MC method, are [λ_1, ..., λ_10, β]=", round.(result_mc, digits=3))

#Calculate estimators, QMC Method
result_qmc=mean.([lambda_qmc[:,j] for j in 1:size(lambda_qmc,2)])
push!(result_qmc, mean(b_qmc))
println("The estimators, calculated using QMC inputs, are [λ_1, ..., λ_10, β]=", round.(result_qmc, digits=3))

#Calculate variances
var_mc=var.([lambda_mc[:,j] for j in 1:size(lambda_mc,2)])
push!(var_mc, var(b_mc))
var_qmc=var.([lambda_qmc[:,j] for j in 1:size(lambda_qmc,2)])
push!(var_qmc, var(b_qmc))



#Create row descriptions (for table)
rows=[]
for i in 1:10 push!(rows, "λ_"*string(i)) end
push!(rows, "β")

#Create header
header=["Pump" "MC" "QMC" "Ratio"]

data=hcat(rows,round.(var_mc, digits=3), round.(var_qmc, digits=3), round.(var_mc./var_qmc, digits=3))
println("Variance of Simulations, m=", m)
pretty_table(data, header)
