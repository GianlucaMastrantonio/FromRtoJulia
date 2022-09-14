module EsempioJulia

##### Packages
# using Pkg
# Pkg.instantiate()
# Pkg.activate("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/julia/EsempioJulia")

using Distributions, Random
using LinearAlgebra, PDMats, StatsBase
using ToggleableAsserts
using ProgressMeter
# using Optim
# using Clustering
# using Parameters
# using SpecialFunctions
# using StatsFuns

# using ThreadsX
# using Impute
#using SpecialFunctions

##### Abstract types

abstract type Sigma2Update end

struct Sigma2Gibbs <: Sigma2Update
    function Sigma2Gibbs()
       
        new()
    end
end
struct Sigma2Metropolis <: Sigma2Update
    function Sigma2Metropolis()
       
        new()
    end
end
struct Sigma2Hamil <: Sigma2Update
    function Sigma2Hamil()
       
        new()
    end
end
struct Sigma2NoUp <: Sigma2Update
    function Sigma2NoUp()
       
        new()
    end
end

##### Include
include(joinpath("model.jl"))
include(joinpath("density.jl"))
include(joinpath("simulations.jl"))
include(joinpath("mcmcsteps/samplemu.jl"))
include(joinpath("mcmcsteps/samplesigma2.jl"))


##### Functions
export    
    sim,
    MCMC



end # module
