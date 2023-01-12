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
    M = normalized_legendre_matrix(t,n; tmin=minimum(t),tmax=maximum(t))

Calculates a matrix with the `n`-th order Legendre coefficients on time points `t`, a vector or range.
The earliest and latest (minimum and maximum) time points will be taken from the edges of `t`.
You can provide these time points, `tmin` and `tmax`, as options.

```juliadoctests
# up to 3rd order coefficients for 5 to 12
julia> M = normalized_legendre_matrix(5:12,3)

# same expression
julia> M = normalized_legendre_matrix([5,6,7,8,9,10,11,12],3)

# as above with time points ranging 1 and 15
julia> M = normalized_legendre_matrix(5:12,3, tmin=1, tmax=15)
```
"""
function normalized_legendre_matrix(t::Union{UnitRange,StepRange,Vector},n; tmin=minimum(t),tmax=maximum(t))
   if !issorted(t)
      throw(ArgumentError("time points not sorted"))
   end
   if tmin > tmax
      throw(ArgumentError("tmin greater than tmax"))
   end
   x = -1 .+ 2*( (t .- tmin)/(tmax-tmin) )
   M = [normalized_legendre(y,k) for y=x,k=0:n]
   return M
end

"""
    vare = hrv_class(ve,rangeset)

Generate a vector containing residual variances from a vector of heterogeneous residual variances `ve` and a set of ranges `rangeset`.

```juliadoctests
julia> vare = hrv_class([10.0,20.0],[3:5,6:9])
7-element Vector{Float64}:
10.0
10.0
10.0
20.0
20.0
20.0
20.0
```
"""
function hrv_class(ve,rangeset)
   if length(ve)!=length(rangeset)
      throw(DimensionMismatch("ve and rangeset"))
   end
   tmin = minimum(rangeset[1])
   tmax = maximum(rangeset[end])
   vare = zeros(tmax-tmin+1)
   for i in 1:length(rangeset)
      for t in rangeset[i]
         vare[t-tmin+1] = ve[i]
      end
   end
   return vare
end
