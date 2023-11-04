
dir_project = "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/gitrepo/FromRtoJulia"
using Pkg
Pkg.activate(dir_project)

using Revise
using Distributions, Random
using LinearAlgebra, PDMats, StatsBase
using RtoJulia
using RCall

#### #### #### #### #### #### 
#### Data simulation
#### #### #### #### #### #### 

rngseed = 123456;
Random.seed!(rngseed);

n::Int64 = 10000
d::Int64 = 2

xmat = Matrix{Float64}(undef, n,d)
for i = 1:d
    xmat[:,i] = rand(Normal(0.0,2.0), n)
end

@rput xmat;

R"""
plot(xmat[,1])
"""

sigma2::Float64 = 3.0
regcoef::Vector{Float64} = rand(Normal(0.0,1.0),d);
ysim::Vector{Float64} = zeros(Float64,n)
xbeta = xmat*regcoef
for i = 1:n
    ysim[i] = rand(Normal(xbeta[i], sigma2^0.5))
end



regcoefOUT, sigma2OUT = mcmc(
    ysim,
    xmat,
    (iter = 10000, thin = 2, burnin = 2),
    (regcoef=[0.0,1000.0], sigma2 = [1.0,1.0]),
    (regcoef = [0.0 for i = 1:size(xmat,2)], sigma2 = [1.0]),
)

@rput regcoefOUT;
@rput sigma2OUT;
@rput regcoef;
@rput sigma2;

R"""
par(mfrow=c(3,3))
i =1
plot(regcoefOUT[i,], type="l")
abline(h =regcoefOUT[i], col=2, lwd=2)
acf(regcoefOUT[i,])
hist(regcoefOUT[i,])

i =2
plot(regcoefOUT[i,], type="l")
abline(h =regcoefOUT[i], col=2, lwd=2)
acf(regcoefOUT[i,])
hist(regcoefOUT[i,])


plot(sigma2OUT[1,], type="l")
abline(h =sigma2, col=2, lwd=2)
acf(sigma2OUT[1,])
hist(sigma2OUT[1,])

"""










