using CSV, Distributions, LinearAlgebra, PrettyTables

#Import the data  (https://fred.stlouisfed.org/graph/?g=r6hR)
Y=CSV.read("/Users/JacobRaymond 1/Library/Mobile Documents/com~apple~CloudDocs/Maitrise/Papier/MCQMC/WTB6MS.csv").WTB6MS./100
Y=Y[497:520]
Y=Y./12

r_iterations=[]
Sige_iterations=[]
aq_iterations=[]
bq_iterations=[]
bp_iterations=[]
ap_iterations=[]
sig_iterations=[]

#Set variables
t=length(Y)
Sig=[100 0.0; 0.0 100]
Sig=inv(Sig)
mu_p=[0.00, 0.00]

#Initial values
r=Y
sig=rand(InverseGamma(1))
a_q=0.01/12
b_q=0.01/12
a_p=0.0001/12
b_p=0.0001/12

for j in 1:300

    #Vector to house the observations
    ap_gibbs=[]
    bp_gibbs=[]
    aq_gibbs=[a_q]
    bq_gibbs=[b_q]
    sig_gibbs=[sig]
    Sige_gibbs=[]
    r_gibbs=[r]


    #Beta distributions
    beta_y=function(a, b, sig, tau)
        gam=sqrt(b^2+2*sig)
        (a/sig)*(2*log(2*gam/((b+gam)*(exp(gam*tau)-1)+2*gam))+(b+gam)*tau)
    end

    beta_r=function(b, sig, tau)
        gam=sqrt(b^2+2*sig)
        2*(1-exp(gam*tau))/((b+gam)*(exp(gam*tau)-1)+2*gam)
    end

    for i in 1:1000
        global sig
        global r
        global b_q
        global a_q
        global b_p
        global a_p

        #Create matrix components for r
        l11=sum(1 ./(sig.*r))
        l12= (t-1)/sig
        l22= sum(r./sig)
        A=sum(diff(r)./r[1:length(r)-1])
        B=sum(-diff(r))

        #Find the mean
        mu=[(l22*A-l12*B)/(l11*l22-l12^2), (l11*B-l12*A)/(l11*l22-l12^2)]

        #Variance-covariance matrix
        vcov= [l11 l12; l12 l22]
        vcov=inv(vcov)

        #New (a,b)^P
        mean_p=inv(Sig+vcov)*(vcov*mu+Sig*mu_p)
        abp=rand(MvNormal(mean_p, vcov+Sig))
        if sum((diff(r).-abp[1].-abp[2].*r[1:t-1])./r[1:t-1])>0
            a_p=abp[1]
            b_p=abp[2]
        end
        push!(ap_gibbs, a_p)
        push!(bp_gibbs, b_p)

        #Generate Sig_e
        w=Y.-beta_y(a_q, b_q, sig, t).-beta_r(b_q, sig, t).*r
        Sig_e=rand(InverseWishart(t+1.0, w*transpose(w)+Diagonal(ones(t))))
        push!(Sige_gibbs, Sig_e)

        #Generate the (a,b)^Q
        ab_can=rand(MvNormal([a_q, b_q], [0.001 0.0; 0.0 0.001]))
        a_can=ab_can[1]
        b_can=ab_can[2]
        if all(i -> i > 0, ab_can)
            q_ratio=pdf(MvNormal(beta_y(a_can, b_can, sig, t).-beta_r(b_can, sig, t).*r,Sig_e),Y)/pdf(MvNormal(beta_y(a_q, b_q, sig, t).-beta_r(b_q, sig, t).*r,Sig_e),Y)
            if q_ratio*(pdf(MvNormal([a_q, b_q], [0.001 0.0; 0.0 0.001]),[a_q, b_q])/pdf(MvNormal([a_q, b_q], [0.001 0.0; 0.0 0.001]),ab_can))>rand()
                a_q=a_can
                b_q=b_can
            end
        end
        push!(aq_gibbs, a_q)
        push!(bq_gibbs, b_q)


        #Generate new sigma
        phi_sig=(diff(r).-a_p.-b_p.*r[1:t-1])./r[1:t-1]
        sig_can=rand(InverseGamma(1.0+t/2, 1.0+sum(phi_sig)/2))
        sig_ratio=pdf(MvNormal(beta_y(a_q, b_q, sig_can, t).-beta_r(b_q, sig_can, t).*r,Sig_e),Y)/pdf(MvNormal(beta_y(a_q, b_q, sig, t).-beta_r(b_q, sig, t).*r,Sig_e),Y)
        if rand()< sig_ratio*pdf(InverseGamma(1.0+t/2, 1.0+sum(phi_sig)/2), sig)/pdf(InverseGamma(1.0+t/2, 1.0+sum(phi_sig)/2), sig_can)
            sig=sig_can
        end

        #Generate new r
        r_can=rand(MvNormal(beta_y(a_q, b_q, sig, t).-beta_r(b_q, sig, t).*r,Sig_e))
        if all(i -> i > 0, r_can)
            pr_can=(1/sqrt(prod(r_can)))*exp(-0.5*sum((diff(r).-a_p.-b_p.*r_can[1:t-1])./r_can[1:t-1])/sig)
            pr=(1/sqrt(prod(r)))*exp(-0.5*sum((diff(r).-a_p.-b_p.*r[1:t-1])./r[1:t-1])/sig)
            r_ratio=(pr_can/pr)*pdf(MvNormal(beta_y(a_q, b_q, sig, t).-beta_r(b_q, sig, t).*r,Sig_e), Y)/pdf(MvNormal(beta_y(a_q, b_q, sig, t).-beta_r(b_q, sig, t).*r_can,Sig_e), Y)
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

header=["r", "∑_e", "a^P", "b^P", "a^Q", "b^Q", "sigma^2"]
data=hcat(mean(var(r_iterations)),mean(var(Sige_iterations)), var(ap_iterations), var(bp_iterations), var(aq_iterations), var(bq_iterations), var(sig_iterations))
pretty_table(data, header)
