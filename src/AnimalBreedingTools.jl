__precompile__()

module AnimalBreedingTools

using Printf
using Statistics
using LinearAlgebra
using Base.Threads
using StatsBase
using Jacobi

struct AnimalBreedingDataSetRaw
   pedlist::Matrix{Int}
   genid::Vector{Int}
   G::Matrix{Float64}
   y::Vector{Float64}
   X::Matrix{Float64}
end

export get_design_matrix, directsum, solve_pcg, approximate_trace_of_inverse
export v2r, c2r, r2v, r2c, vech, bending, bending2
export AnimalBreedingDataSetRaw, load_data_set
export chol_symtrid, ldlt_symtrid, takahashi_ldlt_symtrid!
export normalized_legendre, normalized_legendre_matrix, hrv_class
export disp

include("designmatrix.jl")
include("solver.jl")
include("tools.jl")
include("data.jl")
include("polynomials.jl")
include("io.jl")

end
