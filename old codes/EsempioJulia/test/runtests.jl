using EsempioJulia
using Test


function ff1(app::Vector{Float64})

app[1] = 2.2

end

function ff2(app::Vector{Float64})

    app[2] = 2.2
    
end

function ff1(app::Vector{Float64}, mat::Vector{Float64})

    return app === mat
    
end
    
function ff2(app::Vector{Float64},mat::Vector{Float64})

    return app === mat
    
end

Mat = zeros(Float64,10)

ff1(Mat[:]);
ff2(reshape(Mat, (10)));

@testset "Test Base" begin
    @test π ≈ 3.14 atol=0.01
    @test Mat[1] == 2.2
    @test Mat[2] == 2.2
    @test ff1(Mat[:], reshape(Mat, (10)));
end;
