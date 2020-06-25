using LinearAlgebra, Distributions

#u must be of length d

function MVN_QMC(u, mu, Sigma)

    #Cholesky decomposition
    A=cholesky(Sigma)
    A=A.U

    #Generate standard normals using qmc
    z=map(x->quantile(Normal(), x), u)

    #Output result
    mu+A*z
end
