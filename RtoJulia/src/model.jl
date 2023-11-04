
function mcmc(
    data::Vector{Float64},
    X::Matrix{Float64},
    mcmc::NamedTuple{(:iter, :thin, :burnin), Tuple{Int64, Int64, Int64}},
    priors::NamedTuple{(:regcoef, :sigma2), Tuple{Vector{Float64}, Vector{Float64}}},
    init::NamedTuple{(:regcoef, :sigma2), Tuple{Vector{Float64}, Vector{Float64}}}
)::Tuple{Matrix{Float64}, Matrix{Float64}}

    ncov = size(X,2)
    # MCMC Parameters
    
    SampleToSave::Int64 = Int64(trunc(mcmc.iter-mcmc.burnin)/mcmc.thin)
    
    regcoefMCMC = Matrix{Float64}(undef,  ncov,1)
    regcoefMCMC[:] =  init.regcoef
    regcoefOUT = Matrix{Float64}(undef,  ncov, SampleToSave)
    
    sigma2MCMC = Matrix{Float64}(undef,  1,1)
    sigma2MCMC[:] =  init.sigma2
    sigma2OUT = Matrix{Float64}(undef,  1, SampleToSave)
    
    ### MCMC
    thinburnin = mcmc.burnin
      
    for iMCMC = 1:SampleToSave
    
        for _ = 1:thinburnin
                        
            ### sample Mu
            sample_beta!(data,X,  priors.regcoef,regcoefMCMC, sigma2MCMC)
            
            ### sample sigma2
            sample_sigma2!(data,X,  priors.sigma2,regcoefMCMC, sigma2MCMC)
            
            

            
        end
        thinburnin = mcmc.thin
        
        regcoefOUT[:,iMCMC] .= regcoefMCMC
        sigma2OUT[:,iMCMC] .= sigma2MCMC
    end


    return regcoefOUT, sigma2OUT

end