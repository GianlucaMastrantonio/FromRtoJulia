
function MCMC(
    data::Matrix{Float64};
    X::Matrix{Float64},
    mcmc::NamedTuple{(:iter, :thin, :burnin), Tuple{Int64, Int64, Int64}} = (iter = 10000, thin = 2, burnin = 2),
    priors::NamedTuple{(:regcoef, :sigma2), Tuple{Vector{Float64}, Vector{Float64}}} = (regcoef=[0.0,1000.0], sigma2 = [0.1,0.1]),
    init::NamedTuple{(:regcoef, :sigma2), Tuple{Vector{Float64}, Vector{Float64}}},
    sigma2Update::Sigma2Update = Sigma2Gibbs(),
    sigma2distr::Distribution{Univariate, Continuous} = InverseGamma(),
    var_prop::Float64 = 0.1,
    

)::Tuple{Matrix{Float64}, Matrix{Float64}}

    n = size(data,1)
    d = size(data,2)
    ntot::Int64 = n*d
    
    ncov = size(X,2)
    
    # controllo init
    @toggled_assert size(init.regcoef,1) == ncov "beta - numero di coefficienti sbagliati"
    @toggled_assert size(init.sigma2,1) == 1    "sigma2 - numero di coefficienti sbagliati"

    
    # MCMC Parameters
    
    SampleToSave::Int64 = Int64(trunc(mcmc.iter-mcmc.burnin)/mcmc.thin)
    
    regcoefMCMC = Matrix{Float64}(undef,  ncov,1)
    regcoefMCMC[:] =  init.regcoef
    regcoefOUT = Matrix{Float64}(undef,  ncov, SampleToSave)
    
    sigma2MCMC = Matrix{Float64}(undef,  1,1)
    sigma2MCMC[:] =  init.sigma2
    sigma2OUT = Matrix{Float64}(undef,  1, SampleToSave)
    

    ### prior sigma

    if sigma2distr == InverseGamma()
        sigma2distr_updated = InverseGamma(priors.sigma2...)
    end
    
    
    if sigma2distr == Normal()
        # this is actually a truncated normal prior
        sigma2distr_updated = Normal(priors.sigma2...)
    end
    
    if sigma2distr == Gamma()
        sigma2distr_updated = Gamma(priors.sigma2...)
    end
    

    ### MCMC

    thinburnin = mcmc.burnin
    iter = Int64(0)
    p2 = Progress(mcmc.burnin +(SampleToSave-1)*mcmc.thin, desc="iterations ", offset=0,showspeed=true)
        
    ##### 
    
    println("Iterations: ",mcmc.iter)
    println("burnin: ",mcmc.burnin)
    println("thin: ",mcmc.thin)
    println("number of posterior samples: ",SampleToSave)
    println("Number of threads: ",Threads.nthreads())

    for iMCMC = 1:SampleToSave
    
        for _ = 1:thinburnin
        
            iter += 1
            ProgressMeter.next!(p2; showvalues = [(:iterationstot,mcmc.iter), (:iterations,iter)])
            
            ### sample Mu
            samplemu!(data,X,  priors.regcoef,regcoefMCMC, sigma2MCMC)
            
            ### sample sigma2
            samplesigma2!(data,X,  priors.sigma2,regcoefMCMC, sigma2MCMC, sigma2Update, var_prop)
            #samplesigma2!(data,X,  priors.sigma2,regcoefMCMC, sigma2MCMC, sigma2distr_updated, var_prop)
            
            ### adaptive sigma2
            adaptive_sigma(sigma2Update, iter)
            
        end
        thinburnin = mcmc.thin
        
        regcoefOUT[:,iMCMC] .= regcoefMCMC
        sigma2OUT[:,iMCMC] .= sigma2MCMC
    end


    return regcoefOUT, sigma2OUT

end