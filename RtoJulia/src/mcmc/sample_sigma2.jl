
function sample_sigma2!(data::Vector{Float64},X::Matrix{Float64},  prior::Vector{Float64},regcoefMCMC::Matrix{Float64}, sigma2MCMC::Matrix{Float64})

    ncov = size(X,2)

    mean = X*regcoefMCMC[:]
    
    para::Float64 = Float64(size(data[:],1)/2.0) +  prior[1]
    parb::Float64 = transpose((data[:]-mean))*(data[:]-mean)/2.0 +  prior[2]
    
    sigma2MCMC[:] .= 1/rand(Gamma(para,1.0/parb))
    

end
