module RtoJulia

    using Distributions
    using Random
    using LinearAlgebra
    using PDMats
    using StatsBase

    include(joinpath("model.jl"))
    include(joinpath("mcmc/sample_beta.jl"))
    include(joinpath("mcmc/sample_sigma2.jl"))

    export mcmc

end # module
