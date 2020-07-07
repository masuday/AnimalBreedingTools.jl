__precompile__()

module AnimalBreedingTools

using Printf
using LinearAlgebra
using StatsBase
using Base.Threads

struct AnimalBreedingDataSetRaw
   pedlist::Matrix{Int}
   genid::Vector{Int}
   G::Matrix{Float64}
   y::Vector{Float64}
   X::Matrix{Float64}
end

export get_design_matrix, directsum, solve_pcg, approximate_trace_of_inverse
export v2r, c2r, r2v, r2c, vech
export AnimalBreedingDataSetRaw, load_data_set

include("designmatrix.jl")
include("solver.jl")
include("tools.jl")
include("data.jl")

end
