"""
    v = normalized_legendre(x,n)

Calculates a single value of the `n`-th order Legendre polynomial of `x`.
It calls `Jacobi.legendre`.
"""
function normalized_legendre(x,n)
   a = sqrt((2*n+1)/2)
   return a*legendre(x,n)
end

"""
    M = normalized_legendre_matrix(t,n)

Calculates a matrix with the `n`-th order Legendre coefficients on time points `t`, a vector or range.

```juliadoctests
# up to 3rd order coefficients for 5 to 12
M = normalized_legendre_matrix(5:12,3)

# same expression
M = normalized_legendre_matrix([5,6,7,8,9,10,11,12],3)
```
"""
function normalized_legendre_matrix(t::Union{UnitRange,Vector},n)
   if !issorted(t)
      throw(ArgumentError("time points not sorted"))
   end
   tmin = minimum(t)
   tmax = maximum(t)
   x = -1 .+ 2*( (t .- tmin)/(tmax-tmin) )
   M = [normalized_legendre(y,k) for y=x,k=0:n]
   return M
end
