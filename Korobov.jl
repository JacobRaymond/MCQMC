function Korobov(N, a, d)
    #Create the vector that will house the powers of a
    g=[1]
    for i in 2:d
        push!(g,(a*g[i-1])%N)
    end

    #Initial point of the sequence
    u=[zeros(d)]

    #Korobov sequence
    for j in 1:N-1
        ind=collect((j-1):(j+d-2)).%a .+1
        push!(u, ind.*g./N.%1)
    end

    #Cranley-Patterson Rotation
    ucp=[rand(d)]
    map.([x->x%1], (u.+ucp))
end
