
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

"""
   A = unvech(v; n)

Convert a half-stored vector `v` to a square matrix `A` with order `n`.
It is the inverse function of `vech`.
When `n` is missing, the function will guess the order.
"""
function unvech(v, n=nothing)
   if(isnothing(n))
      a = length(v)
      m = Int.((-1+sqrt(1+8*a))/2)
   else
      m = n
   end
   M = zeros(m,m)
   k = 0
   for j in 1:m
      for i in j:m
         k = k + 1
         M[i,j] = v[k]
         M[j,i] = v[k]
      end
   end
   return(M)
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

"""
    bentB = bending(B [, W], gamma=nothing, tol=eps, verbose=false)
    bentB = bending(B, gamma=nothing, tol=eps, verbose=false)

Apply a traditional "bending" method to a symmetric matrix `B` (between-group-variance) using another symmetric matrix `W` (within-group-variance), as described by Hayes and Hill (1981).
The resulting matrix is `bentB = (1-gamma)*B + gamma*v*W`, where `v` is the average eigenvalues of `B` and `gamma` is a constant between 0 and 1, which can be provided by the user or generated by this function.
When `W` is omitted, the identity matrix will be used instead.

If the user supplies `gamma`, the function will simply apply the bending formula and return `bentB` without checking whether `bentB` is positive defnite or not.
Otherwise (`gamma=nothing`), this function generates `gamma=1/1000`, applies bending, and returns `bentB` if it is positive definite as `minimum(eigen(bentB).values)>tol` (the default `tol` is the machine epsilon).
If not, the function tries `gamma=2/1000`, `gamma=3/1000`, ..., and `gamma=1000/1000` until `bentB` becomes positive definite.

With `verbose=true`, this function prints `v` (average eigenvalue of `B`) and `gamma` that has been used for bending (default: `verbose=false`).
"""
function bending(B::Matrix{Tv}, W=I(size(B,1)); gamma::Union{Nothing,Float64}=nothing, verbose::Bool=false, tol=eps(one(Tv))) where Tv<:AbstractFloat
   if !issymmetric(B)
      throw(ArgumentError("B must be symmetric."))
   end
   if !issymmetric(W)
      throw(ArgumentError("W must be symmetric."))
   end

   v = mean(eigen(B).values)
   if verbose; println("mean(eig(B)): $(v)"); end
   if !isnothing(gamma)
      if gamma<0 || gamma>1
         throw(ArgumentError("gamma should be between 0 and 1."))
      else
         if verbose; println("gamma: $(gamma)"); end
         M = (1-gamma)*B + gamma*v*W
         return M
      end
   end

   M = copy(B)
   if minimum(eigen(M).values)>tol
      return M
   end

   for i=1000:-1:1
      g = 1/i
      M .= (1-g)*B + g*v*W
      if minimum(eigen(M).values)>tol
         if verbose; println("gamma: $(g)"); end
         return M
      end
   end

   return M
end

"""
    bentV = bending2(V [, W], corr=false, mineig=eps, tol=eps, maxiter=10000, force=false, verbose=false)

Apply an alternative "bending" method to a variance-covariance matrix `V` using a weighting matrix `W`, as described by Jorjani et al. (2003).
When omitting `W`, `J` (matrix of ones) will be used, implying the same weights for all elements.
Jorjani et al. (2003) suggested `W` as the reciprocal of number of individuals used in the estimation of covariances.

The function calculates the eigenvalues of `V`, and an eigenvalue is replaced with `mineig` if it is smaller than `mineig`.
The covariance matrix is reconstructed with modified eigenvalues, then tests whether it is positive definite by checking the minimum eigenvalue larger than `tol`.
The function repeats the process until the resulting matrix becomes positive definite, and the number of maximum iterations is defind by `maxiter`.

See some options below.

- With `corr=true`, this function assumes the input is a correlation matrix; otherwise, the input is expected to be a variance-covariance matrix.
- With `force=true`, this function iterates the formula `maxiter` times regardless of whether it has converged or not.
- With `verbose=true`, this function shows the number of iterations required when `force=false`.
"""
function bending2(V::Matrix{Tv}, W=ones(size(V,1),size(V,2)); mineig=eps(one(Tv)), tol=eps(one(Tv)), corr=false, maxiter=10000, force=false, verbose=false) where Tv<:AbstractFloat
   if !issymmetric(V)
      throw(ArgumentError("B must be symmetric."))
   end
   if !issymmetric(W)
      throw(ArgumentError("W must be symmetric."))
   end
   if size(V) != size(W)
      throw(DimensionMismatch("V and W"))
   end

   if corr
      dim = size(V,1)
      Rprev = copy(V)
      Rn = copy(V)
      Rnew = copy(V)
      Dn = zeros(Tv,dim)
      Δn = zeros(Tv,dim)
      for i=1:maxiter
         (D, U) = eigen(Rprev)
         Dn .= D
         Dn[ D .< tol ] .= 2*mineig
         Δn .= Dn * (sum(D)/sum(Dn))
         Rn .= Rprev - ((Rprev - U*diagm(Δn)*U') .* W)
         Rnew = v2r(Rn)
         if !force && minimum(eigen(Rnew).values)>tol
            if verbose; println("round: $(i)"); end
            break
         end
         Rprev .= Rnew
      end
      return Rnew
   else
      Vprev = copy(V)
      Vnew = copy(V)
      for i=1:maxiter
         (D, U) = eigen(Vprev)
         D[ D .< tol ] .= mineig
         Vnew .= Vprev - ((Vprev - U*diagm(D)*U') .* W)
         if !force && minimum(eigen(Vnew).values)>tol
            if verbose; println("round: $(i)"); end
            break
         end
         Vprev .= Vnew
      end
      return Vnew
   end
end
