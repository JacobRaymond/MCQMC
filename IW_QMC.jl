using LinearAlgebra, Distributions

#u must be of length ∑^d_{i=1} i

IW_QMC=function(u, X, v)
    #Dimensions
    d=size(X,1)

    #Create Diagonal
    c=Float64[]
    for i in 1:d
        push!(c, sqrt(quantile(Chisq(v+1-i), u[i])))
    end

    #"A" matrix
    n=map(x->quantile(Normal(), x), u[d+1:length(u)])
    A=[x>y ? pop!(n) : 0 for x in 1:d, y in 1:d]
    A=A+Diagonal(c)

    #Bartlett Decomposition
    D=cholesky(Hermitian(inv(X))).L
    DA=inv(D*A)
    Transpose(DA)*DA
end
