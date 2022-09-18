

function densy(data::Vector{Float64}, mean::Vector{Float64}, var::Float64)::Float64

    n::Int64 = size(data,1)
    ret = zeros(Float64, n)
    for iobs = 1:n
        ret[iobs] = logpdf(Normal(mean[iobs], var^0.5),data[iobs])
    end
    
    # Threads.@threads for iobs = 1:n
    #     println(Threads.threadid())
    #     ret[iobs] = logpdf(Normal(mean[iobs], var^0.5),data[iobs])
    # end
    
    # @sync  for iobs = 1:n
    #     Threads.@spawn ret[iobs] = logpdf(Normal(mean[iobs], var^0.5),data[iobs])
    # end
    
    # @sync begin
    #     for iobs = 1:n
    #         Threads.@spawn begin
    #             println(Threads.threadid())
    #             ret[iobs] = logpdf(Normal(mean[iobs], var^0.5),data[iobs])
    #         end
             
    #     end
    # end
    
    return sum(ret)

end



function densy_spawn(data::Vector{Float64}, mean::Vector{Float64}, var::Float64)

    n::Int64 = size(data,1)
    ret = zeros(Float64, n)
    id = zeros(Int64, n,2)
    @sync begin
        for iobs = 1:n
            Threads.@spawn begin
                id[iobs, :] = [iobs,Threads.threadid() ]
                ret[iobs] = logpdf(Normal(mean[iobs], var^0.5),data[iobs])
                            end
            
        end
    end
    
    return id

end

function densy_thread(data::Vector{Float64}, mean::Vector{Float64}, var::Float64)

    n::Int64 = size(data,1)
    ret = zeros(Float64, n)
    id = zeros(Int64, n,2)
    
    
    Threads.@threads for iobs = 1:n
        id[iobs, :] = [iobs,Threads.threadid() ]
        ret[iobs] = logpdf(Normal(mean[iobs], var^0.5),data[iobs])
    end
    
    
    return id

end

