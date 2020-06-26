__precompile__()

module AnimalBreedingTools

using Printf
using LinearAlgebra

export get_design_matrix, directsum, solve_pcg

include("designmatrix.jl")
include("solver.jl")

end
