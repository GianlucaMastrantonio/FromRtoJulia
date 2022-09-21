

function adaptive_sigma(sigma2Update::Sigma2Update, iter::Int64)


end

function adaptive_sigma(sigma2Update::Sigma2AdaptMetropolis, iter::Int64)

    if mod(iter,sigma2Update.iterupdate) == 0
        
        
        
        mean_alpha = sigma2Update.sumalpha[1]/sigma2Update.iterupdate 
        
        #println(sigma2Update.var[1], " ", mean_alpha)
    
        sigma2Update.var[1] = exp(log(sigma2Update.var[1])+ 100.0/(200.0+iter)*(mean_alpha - sigma2Update.alphatarget)   )
    
        sigma2Update.sumalpha[1] = 0.0
        
        #println(sigma2Update.var[1])
    
        
    
    end

end



function samplesigma2!(data::Matrix{Float64},X::Matrix{Float64},  prior::Vector{Float64},regcoefMCMC::Matrix{Float64}, sigma2MCMC::Matrix{Float64},
    sigma2Update::Sigma2AdaptMetropolis, var_prop::Float64 )

    
    
    prop = rand(Normal(log(sigma2MCMC[1]), sigma2Update.var[1]^0.5))
    var_prop = exp(prop)
    
    Xbeta::Vector{Float64} = X*regcoefMCMC[:]
    
    logdens_prop = densy(data[:], Xbeta[:], var_prop)
    logdens_acc = densy(data[:], Xbeta[:], sigma2MCMC[1])
    
    prior_prop = -(prior[1]+1)*log(var_prop) - prior[2]/var_prop +log(var_prop)
    prior_acc = -(prior[1]+1)*log(sigma2MCMC[1]) - prior[2]/sigma2MCMC[1] +log(sigma2MCMC[1])
    
    MHalpha = exp(logdens_prop+prior_prop-logdens_acc-prior_acc)
    
    if rand(Uniform(0.0,1.0))< MHalpha
    
        sigma2MCMC[1] = var_prop
    
    end
    
    sigma2Update.sumalpha[1] = sigma2Update.sumalpha[1]+min(1.0,MHalpha)
    
    return nothing

end

function samplesigma2!(data::Matrix{Float64},X::Matrix{Float64},  prior::Vector{Float64},regcoefMCMC::Matrix{Float64}, sigma2MCMC::Matrix{Float64},
    sigma2Update::Sigma2Metropolis, var_prop::Float64 )

    
    
    prop = rand(Normal(log(sigma2MCMC[1]), var_prop^0.5))
    var_prop = exp(prop)
    
    Xbeta::Vector{Float64} = X*regcoefMCMC[:]
    
    logdens_prop = densy(data[:], Xbeta[:], var_prop)
    logdens_acc = densy(data[:], Xbeta[:], sigma2MCMC[1])
    
    prior_prop = -(prior[1]+1)*log(var_prop) - prior[2]/var_prop +log(var_prop)
    prior_acc = -(prior[1]+1)*log(sigma2MCMC[1]) - prior[2]/sigma2MCMC[1] +log(sigma2MCMC[1])
    
    MHalpha = exp(logdens_prop+prior_prop-logdens_acc-prior_acc)
    
    if rand(Uniform(0.0,1.0))< MHalpha
    
        sigma2MCMC[1] = var_prop
    
    end
    
    return nothing

end


function samplesigma2!(data::Matrix{Float64},X::Matrix{Float64},  prior::Vector{Float64},regcoefMCMC::Matrix{Float64}, sigma2MCMC::Matrix{Float64},
    sigma2Update::Sigma2NoUp, var_prop::Float64 )

    

end

function samplesigma2!(data::Matrix{Float64},X::Matrix{Float64},  prior::Vector{Float64},regcoefMCMC::Matrix{Float64}, sigma2MCMC::Matrix{Float64},
    sigma2Update::Sigma2Hamil, var_prop::Float64 )

    println("Quando mi va la implemento")
    error("")
    

end

function samplesigma2!(data::Matrix{Float64},X::Matrix{Float64},  prior::Vector{Float64},regcoefMCMC::Matrix{Float64}, sigma2MCMC::Matrix{Float64},
    sigma2Update::Sigma2Gibbs, var_prop::Float64
    )

    ncov = size(X,2)

    mean = X*regcoefMCMC[:]
    
    para::Float64 = Float64(size(data[:],1)/2.0) +  prior[1]
    parb::Float64 = transpose((data[:]-mean))*(data[:]-mean)/2.0 +  prior[2]
    
    sigma2MCMC[:] .= rand(Gamma(para,1.0/parb))
    

end