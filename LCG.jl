using Random

function lcg(N, a, d, seed)


    #Generate starting value
    u=[seed]
    x=N*u

    #Generate the x and u
    for i in 1:N
        x_i=(a*x[i])%N
        push!(x, x_i)
        push!(u, x_i/N)
    end

    #Remove first observations
    popfirst!(x)
    popfirst!(u)

    #Initial point of the sequence
    z=[zeros(d)]

    #Create points z
    for j in 0:N-1
        ind=collect((d*j+1):(d*j+d))
        ind=(ind.-1).%N .+1
        push!(z, u[ind])
     end

    #Cranley-Patterson Rotation
    ucp=rand(d)
    rqmc=map(x->(x.+ucp).%1, z)

end
