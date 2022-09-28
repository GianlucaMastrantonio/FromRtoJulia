# function sima(X, reg,s1)

#     ret = zeros(Float64,size(X,1) )
    
#     Xbeta = X*reg
    
#     for i = axes(X,1)
    
#         ret[i] = rand(Normal(Xbeta[i], s1^0.5  ))
#     end
#     return ret
# end
function sim(X::Matrix{Float64}, reg::Vector{Float64},s1::Float64)::Vector{Float64}

    ret = zeros(Float64,size(X,1) )
    
    Xbeta = X*reg
    
    for i = axes(X,1)
    
        ret[i] = rand(Normal(Xbeta[i], s1^0.5  ))
    end
    return ret
end

function sim(ret::Vector{Float64},X::Matrix{Float64}, reg::Vector{Float64},s1::Float64)

    
    Xbeta = X*reg
    
    for i = axes(X,1)
    
        ret[i] = rand(Normal(Xbeta[i],  s1^0.5 ))
    end
    return nothing
end
