using Random

function lcg(N, a, d)


    #Generate starting value
    u=[1/N]
    x=N*u

    #Generate the x and u
    for i in 1:N-2
        x_i=(a*x[i])%N
        push!(x, x_i)
        push!(u, x_i/N)
    end

    #Initial point of the sequence
    z=[zeros(d)]

    #Generate the order of the first indexes for each point
    seq=[1:1:N-1;]
    indices=[1]
    for j in 1:N-2
        ind=(indices[j]+d)%(N-1)
        if ind==0 || ind in indices
            ind=minimum(setdiff(seq, indices))
        end
        push!(indices, ind)
    end

    #Order points
    indices=map(x->[x:1:x+d-1;], indices)
    indices=map(x->(x.-1).%(N-1).+1, indices)

    #Create Points
    for i in 1:length(indices)
        push!(z, u[indices[i]])
    end

    #Cranley-Patterson Rotation
    ucp=rand(d)
    rqmc=map(x->(x.+ucp).%1, z)

end
