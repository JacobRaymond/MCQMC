using Distributions

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
lambda_gibbs=[ones(10)]
b_gibbs=Float16[]

#Gibbs Sampler

n=300 #Number of times the sampler is repeated
for k in 1:n
    m=1021 #Number of iterations of the sampler
    for j in 1:m
        global b
        global lambda
        global s
        global t

        #Draw the lambdas
        lambda=Float16[rand(Gamma(a+s[1], 1/(t[1]+1/b)))]
        for i in 2:10
            push!(lambda, rand(Gamma(a+s[i], 1/(t[i]+1/b))))
        end

        #Draw the beta
        b=rand(InverseGamma(gam+10*a, sum(lambda)+del))
    end

    #Save the values
    push!(lambda_gibbs, lambda)
    push!(b_gibbs, b)
end
#Calcualte the results
deleteat!(lambda_gibbs, 1) #The first observation of lambda_gibbs is an array of 1s due to storage issues
result=mean(lambda_gibbs)
push!(result, mean(b_gibbs))
println("[λ_1, ..., λ_10, β]=", round.(result, digits=3))

#### RQMC Method (see  DOI: 10.1073/pnas.0409596102)####

#Import Korobov sampler
include("Korobov.jl")

#Initial values
init=s./t #MLE
b=(sum(init)+del)/(gam+10*a-1)#Mean of a IG(gamma+10a, ∑lambda+∂)

#Array to house the observations
lambda_gibbs=[ones(10)]
b_gibbs=Float16[]

#Draw the lambdas
m=200 #The algorithm is repeated 300 times, as per Owen and Tribble
for k in 1:m
    #Generate RMQC points
    rqmc=Korobov(1021, 65, 11)

    for j in 1:1021
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
    push!(lambda_gibbs, lambda)
    push!(b_gibbs, b)
end

#Calcualte the results
deleteat!(lambda_gibbs, 1) #The first observation of lambda_gibbs is an array of 1s due to storage issues
result=mean(lambda_gibbs)
push!(result, mean(b_gibbs))
println("[λ_1, ..., λ_10, β]=", round.(result, digits=3))
