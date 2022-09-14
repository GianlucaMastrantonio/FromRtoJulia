

function densy(data::Vector{Float64}, mean::Vector{Float64}, var::Float64)::Float64

    n::Int64 = size(data,1)
    ret = zeros(Float64, n)
    for iobs = 1:n
        ret[iobs] = logpdf(Normal(mean[iobs], var^0.5),data[iobs])
    end
    
    # Threads.@threads for iobs = 1:n
    #     ret[iobs] = logpdf(Normal(mean[iobs], var^0.5),data[iobs])
    # end
    
    # @sync  for iobs = 1:n
    #     Threads.@spawn ret[iobs] = logpdf(Normal(mean[iobs], var^0.5),data[iobs])
    # end
    
    
    return sum(ret)

end

