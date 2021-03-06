using CSV, Distributions, LinearAlgebra, PrettyTables
include("LCG.jl")
include("IW_QMC.jl")
include("MVN_QMC.jl")

#Import the data  (https://fred.stlouisfed.org/graph/?g=r6hR)
Y=CSV.read("/Users/JacobRaymond 1/Library/Mobile Documents/com~apple~CloudDocs/Maitrise/Papier/MCQMC/WTB6MS.csv").WTB6MS./100
Y=Y[310:335]

#### MCMC ####


@time begin

r_iterations=[]
Sige_iterations=[]
aq_iterations=[]
bq_iterations=[]
bp_iterations=[]
ap_iterations=[]
sig_iterations=[]

#Set variables
t=length(Y)
Sig_ab=[0.1 0.0; 0.0 0.1]
Sig=inv(Sig_ab)
mu_p=[0.01, 0.01]
T=26 #Goes into loading functions; 6 months expressed in weeks


for j in 1:100

    #Initial values
    r=Y
    sig=rand()
    a_q=rand()
    b_q=rand()
    a_p=rand()
    b_p=rand()

    #Vector to house the observations
    ap_gibbs=[]
    bp_gibbs=[]
    aq_gibbs=[a_q]
    bq_gibbs=[b_q]
    sig_gibbs=[sig]
    Sige_gibbs=[]
    r_gibbs=[r]


    #Loading functions distributions
    beta_y=function(a, b, sig)
        gam=sqrt(b^2+2*sig)
        (a/sig)*(2*log(2*gam/((b+gam)*(exp(gam*T)-1)+2*gam))+(b+gam)*T)
    end

    beta_r=function(b, sig)
        gam=sqrt(b^2+2*sig)
        2*(1-exp(gam*T))/((b+gam)*(exp(gam*T)-1)+2*gam)
    end

    for i in 1:1021

        #Create matrix components for r
        l11=sum(1 ./(sig.*r[1:length(r)-1]))
        l12= -(t-1)/sig
        l22= sum(r[1:length(r)-1]./sig)
        A=sum(-diff(r)./r[1:length(r)-1])
        B=sum(-diff(r))

        #Find the mean
        mu=[(l22*A-l12*B)/(l11*l22-l12^2), (l11*B-l12*A)/(l11*l22-l12^2)]

        #Variance-covariance matrix
        vcov= [l11 l12; l12 l22]
        vcov=inv(vcov)

        #New (a,b)^P
        Sig_p=inv(Sig+vcov)
        mean_p=Sig_p*(vcov*mu+Sig*mu_p)
        abp=rand(MvNormal(mean_p, Matrix(Hermitian(Sig_p))))
        if sum((diff(r).-abp[1].-abp[2].*r[1:t-1])./r[1:t-1])>0
            a_p=abp[1]
            b_p=abp[2]
        end
        push!(ap_gibbs, a_p)
        push!(bp_gibbs, b_p)

        #Generate Sig_e
        w=Y.-beta_y(a_q, b_q, sig).-beta_r(b_q, sig).*r
        Sig_e=rand(InverseWishart(t+1.0, Matrix(Hermitian(w*transpose(w)+Diagonal(ones(t)))))) #Assuming a prior matrix of I and degs of freedom of t
        push!(Sige_gibbs, Sig_e)

        #Generate the (a,b)^Q
        ab_can=rand(MvNormal([a_q, b_q], [1.0 0.0; 0.0 1.0])) #The identity matrix is relatively fat tailed.
        a_can=ab_can[1]
        b_can=ab_can[2]
        if all(i -> i > 0, ab_can)
            q_ratio=pdf(MvNormal(beta_y(a_can, b_can, sig).-beta_r(b_can, sig).*r,Sig_e),Y)/pdf(MvNormal(beta_y(a_q, b_q, sig).-beta_r(b_q, sig).*r,Sig_e),Y)
            q_ratio=q_ratio*pdf(MvNormal(mu_p, Sig_ab),ab_can)/pdf(MvNormal(mu_p, Sig_ab),[a_q, b_q]) #Prior distribution
            if q_ratio*(pdf(MvNormal([a_q, b_q], [1.0 0.0; 0.0 1.0]),[a_q, b_q])/pdf(MvNormal([a_q, b_q], [1.0 0.0; 0.0 1.0]),ab_can))>rand()
                a_q=a_can
                b_q=b_can
            end
        end
        push!(aq_gibbs, a_q)
        push!(bq_gibbs, b_q)


        #Generate new sigma
        phi_sig=((diff(r).-a_p.-b_p.*r[1:t-1]).^2)./r[1:t-1]
        sig_can=rand(InverseGamma(1.0+(t-1)/2, 1.0+sum(phi_sig)/2))
        sig_ratio=pdf(MvNormal(beta_y(a_q, b_q, sig_can).-beta_r(b_q, sig_can).*r,Sig_e),Y)/pdf(MvNormal(beta_y(a_q, b_q, sig).-beta_r(b_q, sig).*r,Sig_e),Y)
        if rand()< sig_ratio*pdf(InverseGamma(1.0+(t-1)/2, 1.0+sum(phi_sig)/2), sig)/pdf(InverseGamma(1.0+(t-1)/2, 1.0+sum(phi_sig)/2), sig_can)
            sig=sig_can
        end

        #Generate new r
        r_can=rand(MvNormal(r, var(Y).*Matrix(Diagonal(ones(t)))))
        if all(i -> i > 0, r_can)
            pr_can=(1/sqrt(prod(r_can)))*exp(-0.5*sum((diff(r).-a_p.-b_p.*r_can[1:t-1])./r_can[1:t-1])/sig)
            pr=(1/sqrt(prod(r)))*exp(-0.5*sum((diff(r).-a_p.-b_p.*r[1:t-1])./r[1:t-1])/sig)
            r_ratio=(pr_can/pr)*pdf(MvNormal(beta_y(a_q, b_q, sig).-beta_r(b_q, sig).*r_can,Sig_e), Y)/pdf(MvNormal(beta_y(a_q, b_q, sig).-beta_r(b_q, sig).*r,Sig_e), Y)
            r_ratio=r_ratio*(pdf(MvNormal(Y, var(Y).*Matrix(Diagonal(ones(t)))), r))/(pdf(MvNormal(Y, var(Y).*Matrix(Diagonal(ones(t)))), r_can))
            if r_ratio>rand()
                r=r_can
            end
        end
        push!(r_gibbs, r)
    end

    push!(r_iterations, mean(r_gibbs))
    push!(Sige_iterations, mean(Sige_gibbs))
    push!(aq_iterations, mean(aq_gibbs))
    push!(bq_iterations, mean(bq_gibbs))
    push!(ap_iterations, mean(ap_gibbs))
    push!(bp_iterations, mean(bp_gibbs))
    push!(sig_iterations, mean(sig_gibbs))
end
end

#### MCQMC ####

@time begin

r_iterations_qmc=[]
Sige_iterations_qmc=[]
aq_iterations_qmc=[]
bq_iterations_qmc=[]
bp_iterations_qmc=[]
ap_iterations_qmc=[]
sig_iterations_qmc=[]

#Set variables
t=length(Y)
Sig_ab=[0.1 0.0; 0.0 0.1]
Sig=inv(Sig_ab)
mu_p=[0.01, 0.01]
T=26 #Goes into loading functions; 6 months expressed in weeks

for j in 1:100

    #Initial values
    r=Y
    sig=rand()
    a_q=rand()
    b_q=rand()
    a_p=rand()
    b_p=rand()

    #Vector to house the observations
    ap_gibbs=[]
    bp_gibbs=[]
    aq_gibbs=[a_q]
    bq_gibbs=[b_q]
    sig_gibbs=[sig]
    Sige_gibbs=[]
    r_gibbs=[r]

    #Generate qmc points
    u=lcg(1021, 76, Integer(0.5*t^2+1.5*t+8))

    #Loading functions distributions
    beta_y=function(a, b, sig)
        gam=sqrt(b^2+2*sig)
        (a/sig)*(2*log(2*gam/((b+gam)*(exp(gam*T)-1)+2*gam))+(b+gam)*T)
    end

    beta_r=function(b, sig)
        gam=sqrt(b^2+2*sig)
        2*(1-exp(gam*T))/((b+gam)*(exp(gam*T)-1)+2*gam)
    end

    for i in 1:1021

        #Create matrix components for r
        l11=sum(1 ./(sig.*r[1:length(r)-1]))
        l12= -(t-1)/sig
        l22= sum(r[1:length(r)-1]./sig)
        A=sum(-diff(r)./r[1:length(r)-1])
        B=sum(-diff(r))

        #Find the mean
        mu=[(l22*A-l12*B)/(l11*l22-l12^2), (l11*B-l12*A)/(l11*l22-l12^2)]

        #Variance-covariance matrix
        vcov= [l11 l12; l12 l22]
        vcov=inv(vcov)

        #New (a,b)^P
        Sig_p=inv(Sig+vcov)
        mean_p=Sig_p*(vcov*mu+Sig*mu_p)
        abp=MVN_QMC(u[i][1:2], mean_p, Matrix(Hermitian(Sig_p)))
        if sum((diff(r).-abp[1].-abp[2].*r[1:t-1])./r[1:t-1])>0
            a_p=abp[1]
            b_p=abp[2]
        end
        push!(ap_gibbs, a_p)
        push!(bp_gibbs, b_p)

        #Generate Sig_e
        w=Y.-beta_y(a_q, b_q, sig).-beta_r(b_q, sig).*r
        Sig_e=IW_QMC(u[i][3:Integer(0.5*(t^2+t)+2)],w*transpose(w)+Diagonal(ones(t)), t+1.0 ) #Assuming a prior matrix of I and degs of freedom of t
        push!(Sige_gibbs, Sig_e)

        #Generate the (a,b)^Q
        ab_can=MVN_QMC(u[i][Integer(0.5*(t^2+t)+3):Integer(0.5*(t^2+t)+4)],[a_q, b_q], [1.0 0.0; 0.0 1.0]) #The identity matrix is relatively fat tailed.
        a_can=ab_can[1]
        b_can=ab_can[2]
        if all(i -> i > 0, ab_can)
            q_ratio=pdf(MvNormal(beta_y(a_can, b_can, sig).-beta_r(b_can, sig).*r,Sig_e),Y)/pdf(MvNormal(beta_y(a_q, b_q, sig).-beta_r(b_q, sig).*r,Sig_e),Y)
            q_ratio=q_ratio*pdf(MvNormal(mu_p, Sig_ab),ab_can)/pdf(MvNormal(mu_p, Sig_ab),[a_q, b_q]) #Prior distribution
            if q_ratio*(pdf(MvNormal([a_q, b_q], [1.0 0.0; 0.0 1.0]),[a_q, b_q])/pdf(MvNormal([a_q, b_q], [1.0 0.0; 0.0 1.0]),ab_can))>u[i][Integer(0.5*(t^2+t)+5)]
                a_q=a_can
                b_q=b_can
            end
        end
        push!(aq_gibbs, a_q)
        push!(bq_gibbs, b_q)


        #Generate new sigma
        phi_sig=((diff(r).-a_p.-b_p.*r[1:t-1]).^2)./r[1:t-1]
        sig_can=quantile(InverseGamma(1.0+(t-1)/2, 1.0+sum(phi_sig)/2), u[i][Integer(0.5*(t^2+t)+6)])
        sig_ratio=pdf(MvNormal(beta_y(a_q, b_q, sig_can).-beta_r(b_q, sig_can).*r,Sig_e),Y)/pdf(MvNormal(beta_y(a_q, b_q, sig).-beta_r(b_q, sig).*r,Sig_e),Y)
        if u[i][Integer(0.5*(t^2+t)+7)]< sig_ratio*pdf(InverseGamma(1.0+(t-1)/2, 1.0+sum(phi_sig)/2), sig)/pdf(InverseGamma(1.0+(t-1)/2, 1.0+sum(phi_sig)/2), sig_can)
            sig=sig_can
        end

        phi_sig=((diff(r).-a_p.-b_p.*r[1:t-1]).^2)./r[1:t-1]
        sig_can=rand(InverseGamma(1.0+(t-1)/2, 1.0+sum(phi_sig)/2))
        sig_ratio=pdf(MvNormal(beta_y(a_q, b_q, sig_can).-beta_r(b_q, sig_can).*r,Sig_e),Y)/pdf(MvNormal(beta_y(a_q, b_q, sig).-beta_r(b_q, sig).*r,Sig_e),Y)
        if rand()< sig_ratio*pdf(InverseGamma(1.0+(t-1)/2, 1.0+sum(phi_sig)/2), sig)/pdf(InverseGamma(1.0+(t-1)/2, 1.0+sum(phi_sig)/2), sig_can)
            sig=sig_can
        end

        #Generate new r
        r_can=MVN_QMC(u[i][Integer(0.5*(t^2+t)+8):Integer(0.5*(t^2+3*t)+7)],r, var(Y).*Matrix(Diagonal(ones(t))))
        if all(i -> i > 0, r_can)
            pr_can=(1/sqrt(prod(r_can)))*exp(-0.5*sum((diff(r).-a_p.-b_p.*r_can[1:t-1])./r_can[1:t-1])/sig)
            pr=(1/sqrt(prod(r)))*exp(-0.5*sum((diff(r).-a_p.-b_p.*r[1:t-1])./r[1:t-1])/sig)
            r_ratio=(pr_can/pr)*pdf(MvNormal(beta_y(a_q, b_q, sig).-beta_r(b_q, sig).*r_can,Sig_e), Y)/pdf(MvNormal(beta_y(a_q, b_q, sig).-beta_r(b_q, sig).*r,Sig_e), Y)
            r_ratio=r_ratio*(pdf(MvNormal(Y, var(Y).*Matrix(Diagonal(ones(t)))), r))/(pdf(MvNormal(Y, var(Y).*Matrix(Diagonal(ones(t)))), r_can))
            if r_ratio>rand()
                r=r_can
            end
        end
        push!(r_gibbs, r)
    end

    push!(r_iterations_qmc, mean(r_gibbs))
    push!(Sige_iterations_qmc, mean(Sige_gibbs))
    push!(aq_iterations_qmc, mean(aq_gibbs))
    push!(bq_iterations_qmc, mean(bq_gibbs))
    push!(ap_iterations_qmc, mean(ap_gibbs))
    push!(bp_iterations_qmc, mean(bp_gibbs))
    push!(sig_iterations_qmc, mean(sig_gibbs))
end
end

#### Output Results ####

#Show variance for the (a,b)s and sigma^2
header_ab=["","Variance (MC)", "Variance (QMC)", "Ratio"]
cols=["a^P", "b^P", "a^Q", "b^Q", "∑ (Mean)", "σ^2", "r (Mean)"]
data_ab_mc=vcat(var(ap_iterations), var(bp_iterations), var(aq_iterations), var(bq_iterations_qmc), mean(var(Sige_iterations)), var(sig_iterations), mean(var(r_iterations)))
data_ab_qmc=vcat(var(ap_iterations_qmc), var(bp_iterations_qmc), var(aq_iterations_qmc), var(bq_iterations_qmc), mean(var(Sige_iterations_qmc)),var(sig_iterations_qmc), mean(var(r_iterations_qmc)))
data_ab=hcat( data_ab_mc, data_ab_qmc, data_ab_mc./data_ab_qmc)
data_ab=map(x->round.(x, digits=5), data_ab)
data_ab=hcat(cols, data_ab)
pretty_table(data_ab, header_ab)

#Show results for the (a,b)s and sigma^2
header_ab=["","Results (MC)", "Results (QMC)"]
cols=["a^P", "b^P", "a^Q", "b^Q", "σ^2"]
data_ab_mc=vcat(mean(ap_iterations), mean(bp_iterations), mean(aq_iterations), mean(bq_iterations_qmc), mean(sig_iterations))
data_ab_qmc=vcat(mean(ap_iterations_qmc), mean(bp_iterations_qmc), mean(aq_iterations_qmc), mean(bq_iterations_qmc), mean(sig_iterations_qmc))
data_ab=hcat(cols, data_ab_mc, data_ab_qmc)
data_ab[6:15]=map(x->round.(x, digits=5), data_ab[6:15])
pretty_table(data_ab, header_ab)
