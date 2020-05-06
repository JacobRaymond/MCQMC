using Distributions, LinearAlgebra

#### Simple Metropolis-Hastings MCMC ####
function mcmc(niter, mle, prior, init, burnit, param, x, propvar=1)
    acc=0 #Acceptance rate initialization
    est=[init] #Will house the values of the parameter
    propvar=fill(propvar, length(init))
    for i in 1:niter+burnit

        prop=rand(Distributions.MvNormal(init,Diagonal(propvar))) #Generate new proposal value

        if((exp(mle(prop, x)-mle(init, x))+prior(prop)-prior(init))>rand()) #Metropolis algorithm
            init= prop
            if(i>burnit)
                acc +=1
            end
        end
    push!(est, init)
    end

    #Delete burned observations from est
    est=est[(burnit+1):length(est)]

    println(param, " = ", sum(est)/niter)
    println("The acceptance rate is ", 100*acc/niter, "%")
end


#Expand to multiple dimensions


####Simulation Example (One Variable)####
#Simulate some data
var=4
mu=3
n=100
x=(randn(n).*sqrt(var).+mu)

#MLE
function mle1(m, x, v=var)
    (-n/2)*log(v/(2*pi))-sum((x.-m).^2)/(2*v)
end

#Prior (flat)
function prior(y)
    1
end

mcmc(10000, mle1, prior, [0.00], 1000, "mu", x, 1.00)

####Simulation Example (Two Variables)####

#MLE2
function mle2(prop, x)
    if(prop[2]<0) #Account for the possibility that the proposal gives a negative value for sigma
        -Inf
    else
        (-n/2)*log(prop[2]/(2*pi))-sum((x.-prop[1]).^2)/(2*prop[2])
    end
end

mcmc(20000, mle2, prior, [0.00,1.00], 10000, ["mu", "sigma"], x, 0.4) #Larger burn in (50%) because sigma can take negative values.
