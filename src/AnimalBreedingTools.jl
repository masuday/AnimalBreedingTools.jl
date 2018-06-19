__precompile__()

module AnimalBreedingTools

export get_design_matrix, directsum, solve_pcg

include("designmatrix.jl")
include("solver.jl")

end
