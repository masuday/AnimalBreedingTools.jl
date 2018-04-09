using AnimalBreedingTools
using Base.Test

@testset "get_design_matrix" begin
    x = [1,2,1,1,2]
    y = [0.5,0.6,0.7,0.8,0.9]
    # class variable
    X = get_design_matrix(x)
    @test X == [1.0 0.0; 0.0 1.0; 1.0 0.0; 1.0 0.0; 0.0 1.0]
    # covariate
    X = get_design_matrix(x, cov=true)
    @test X == Matrix(convert.(Float64, x[:,:]))
    # nested covariate
    X = get_design_matrix(y, x, cov=true)
    @test X == [0.5 0.0; 0.0 0.6; 0.7 0.0; 0.8 0.0; 0.0 0.9]
end
