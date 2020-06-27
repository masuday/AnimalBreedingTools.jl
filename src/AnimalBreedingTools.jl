__precompile__()

module AnimalBreedingTools

using Printf
using LinearAlgebra

struct AnimalBreedingDataSetRaw
   pedlist::Matrix{Int}
   genid::Vector{Int}
   G::Matrix{Float64}
   y::Vector{Float64}
   X::Matrix{Float64}
end

export get_design_matrix, directsum, solve_pcg
export AnimalBreedingDataSetRaw, load_data_set

include("designmatrix.jl")
include("solver.jl")
include("data.jl")

end
