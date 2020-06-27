
"""
    corr = v2r(cov)

Convert a covariance matrix `cov` to the correlation matrix `corr`.
This is a wrapper of `StatsBase.cov2cor.`
"""
function v2r(C::Matrix{T}) where {T<:AbstractFloat}
   s = sqrt.(diag(C))
   R = StatsBase.cov2cor(C, s)
   return R
end
function c2r(C::Matrix{T}) where {T<:AbstractFloat}
   return v2r(C)
end

"""
    cov = r2v(corr, stdev)

Convert a correlation matrix `cov` and the vector of standard deviation `stdev` to the covariance matrix `cov`.
This is a wrapper of `StatsBase.cor2cov.`
"""
function r2v(C::Matrix{T}, s::Vector{T}) where {T<:AbstractFloat}
   return StatsBase.cor2cov(C,s)
end
function r2c(C::Matrix{T}, s::Vector{T}) where {T<:AbstractFloat}
   return r2v(C,s)
end

"""
   v = vech(A)

Convert a square matrix to a half-stored vector.
It was written by Steven G. Johnson; see `https://discourse.julialang.org/t/half-vectorization/7399` for details.
"""
function vech(A::AbstractMatrix{T}) where T
   m = LinearAlgebra.checksquare(A)
   v = zeros(T,((m*(m+1))>>1))
   k = 0
   for j = 1:m, i = j:m
       @inbounds v[k += 1] = A[i,j]
   end
   return v
end
