# MCQMC
Code to support my master's research on Markov Chain quasi-Monte Carlo methods for term structure model parameter estimation.

## Proofs of Concept
MCMC.jl: A simple Metropolis-Hastings algorithm in Julia, as a proof of concept.

GibbsPoC.jl: An implementation of the second example from "A quasi-Monte Carlo Metropolis algorithm" by Owen and Tribble (DOI: 10.1073/pnas.0409596102). Consists in a Gibbs sampler incorporation RQMC points.

CIRPoC.jl: An implementation of the CIR Gibbs sampler from "Bayesian Estimation of CIR Model" by Feng and Xie (http://www.jds-online.com/files/JDS-746.pdf ) - WORK IN PROGRESS

## QMC Methods
Koborov.jl: Function to generate a Koborov point set.
LCG.jl: Function to generate points from a linear congruential generator.