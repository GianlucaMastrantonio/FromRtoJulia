

function samplemu!(data::Matrix{Float64},X::Matrix{Float64},  prior::Vector{Float64},regcoefMCMC::Matrix{Float64}, sigma2MCMC::Matrix{Float64})

    ncov = size(X,2)

    XtX::Symmetric{Float64, Matrix{Float64}} = Symmetric(transpose(X)*X)
    Xty::Matrix{Float64} = reshape(transpose(X)*data[:], (ncov,1))

    var_app = Symmetric(inv(XtX ./ sigma2MCMC[1] +  Diagonal( [ 1.0/  prior[2] for i = 1:ncov ] )))
    
    
    mean_app = var_app * (Xty ./ sigma2MCMC[1] .+ prior[1] /  prior[2])
    
    #println(var_app[1:3,1:3])
    
    
    @toggled_assert isapprox((inv(XtX ./ sigma2MCMC[1] +  Diagonal( [ 1.0/  prior[2] for i = 1:ncov ] ))),transpose((inv(XtX ./ sigma2MCMC[1] +  Diagonal( [ 1.0/  prior[2] for i = 1:ncov ] )))), rtol=1e-12) "matrice non simmetrica"
    
    regcoefMCMC[:] = rand(MvNormal(mean_app[:],var_app))
    

end