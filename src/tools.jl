
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

# CS 4220/Math 4260 Numerical Analysis: Linear and Nonlinear Problems
# by Charles Van Loan
# http://www.cs.cornell.edu/courses/cs4220/2014sp/CVLBook/chap7.pdf
"""
     dv,ev = chol_symtrid(T::SymTridiagonal)

Calculate the Cholesky factor of a symmetric tridiagonal matrix.
Note that the factor is a lower tridiagonal matrix, and the output is a couple of vectors.
The factor is stored in `dv` for the diagonal matrix (D) and in `ev` for the off-diagonal elements (L).     
"""
function chol_symtrid(T::SymTridiagonal)
   n = length(T.dv)
   dv = zeros(n)
   ev = zeros(n-1)
   dv[1] = sqrt(T.dv[1])
   for i=2:n
      ev[i-1] = T.ev[i-1]/dv[i-1]
      dv[i] = sqrt(T.dv[i] - ev[i-1]^2)
   end
   return dv,ev
end

# The Tridiagonal LDL T Calculus for the Indefinite Generalized Symmetric Eigenproblem
# By Peter Strobach
# https://www.researchgate.net/publication/271196070_THE_TRIDIAGONAL_LDL_T_CALCULUS_FOR_THE_INDEFINITE_GENERALIZED_SYMMETRIC_EIGENPROBLEM
"""
    dv,ev = ldlt_symtrid(T::SymTridiagonal)

Calculate the LDLT factor of a symmetric tridiagonal matrix.
Note that the factor is a lower tridiagonal matrix, and the output is a couple of vectors.
The factor is stored in `dv` for the diagonal matrix (D) and in `ev` for the off-diagonal elements (L).
"""
function ldlt_symtrid(T::SymTridiagonal)
   n = length(T.dv)
   dv = zeros(n)
   ev = zeros(n-1)
   dv[1] = T.dv[1]
   for i=1:n-1
      ev[i] = T.ev[i]/dv[i]
      dv[i+1] = T.dv[i+1] - (ev[i]^2)*dv[i]
   end
   return dv,ev
end

# Takahashi method for the LDLT factor
"""
    takahashi_ldlt_symtrid!(dv,ev)

Compute the sparse inverse (as known as selected inverse) of the LDLT factor of a symmetric trigiagonal matrix.
The LDLT factor is stored in two vectors (`dv` and `ev`), which should be calculated with `ldlt_symtrid`.
"""
function takahashi_ldlt_symtrid!(dv,ev)
   n = size(dv,1)
   dv[n] = 1/dv[n]
   for i=n-1:-1:1
      evi = ev[i]
      ev[i] = -evi*dv[i+1]
      dv[i] = 1/dv[i] - evi*ev[i]
   end
end

function get_invT11_fast(T::SymTridiagonal)
   dv,ev = ldlt_symtrid(T::SymTridiagonal)
   takahashi_ldlt_symtrid!(dv,ev)
   return dv[1]
end