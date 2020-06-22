# MCQMC
Code to support my master's research on Markov Chain quasi-Monte Carlo methods for term structure model parameter estimation.

CIR.jl: A Gibbs sampler to estimate the parameters of a Cox-Ingersoll-Ross model. See http://home.uchicago.edu/~lhansen/JP_handbook.pdf for details - WORK IN PROGRESS

## Proofs of Concept
MCMC.jl: A simple Metropolis-Hastings algorithm in Julia, as a proof of concept.

MCMCPoC.jl: An implementation of the first example from "A quasi-Monte Carlo Metropolis algorithm" by Owen and Tribble (DOI: 10.1073/pnas.0409596102). Consists in a Metropolis-Hastings sampler incorporating RQMC points.

GibbsPoC.jl: An implementation of the second example from "A quasi-Monte Carlo Metropolis algorithm" by Owen and Tribble (DOI: 10.1073/pnas.0409596102). Consists in a Gibbs sampler incorporating RQMC points.

CIRPoC.jl: An implementation of the CIR Gibbs sampler from "Bayesian Estimation of CIR Model" by Feng and Xie (http://www.jds-online.com/files/JDS-746.pdf ). WTB6MS is the required dataset. - WORK IN PROGRESS


MVNPoC.jl: A modification of MCMCPoC.jl in which I check whether Owen and Tribble's method works for multivariate normal distributions.

CIRPoC.R: An implementation of the CIR Gibbs sampler from "Bayesian Estimation of CIR Model" by Feng and Xie (http://www.jds-online.com/files/JDS-746.pdf ). WTB6MS is the required dataset. More advanced than the Julia version, this code iterates the Gibbs sampler multiple times and examines two solutions to prevent negative values for the rates: acceptance/rejection and transformation. R is chosen for two reasons: Efficiency is not emphasized here, and I find it more flexible for experimentation than Julia. - WORK IN PROGRESS

## QMC Methods
Koborov.jl: Function to generate a Koborov point set.

LCG.jl: Function to generate points from a linear congruential generator.

MVN_QMC.jl: Generate QMC draws from a MVN(mu, Sigma) distribution.