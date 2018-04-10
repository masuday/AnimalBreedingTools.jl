"""
    Z = get_design_matrix(x, nested=[]; maxlev=0, cov=false)
    Z = get_design_matrix(X; maxlev=[], cov=[])

Returns a design matrix accoding to class variables or covariates in a vector `x` or
matrix `X`. For a vector `x`, optional argument `dested` defines `x` as the nested
covariate within a class. Option `maxval` specifies the number of levels in this effect.
Omitting `maxval`, the function uses the maximum integer as the maximum level
in `x` or `X`. Option `cov` specifies the input is covariate and passes it through
the output. The matrix `X` doesn't support nested covariates.

```juliadoctest
julia> x = [1,2,1,1,2]
julia> X = get_design_matrix(x)
5×2 Array{Float64,2}:
 1.0  0.0
 0.0  1.0
 1.0  0.0
 1.0  0.0
 0.0  1.0
```
"""
function get_design_matrix(x::Vector{<:Real}, nested::Vector{<:Integer}=Int[]; maxlev=0, cov=false)
   n = length(x)

   # no change
   if length(nested)==0 && cov
      return x[:,:]
   end

   # maximum level
   if maxlev==0
      if length(nested)>0
         lv = maximum(nested)
      else
         lv = maximum(x)
      end
   else
      lv = maxlev
   end
   # fill the value
   Z = zeros(Float64,n,lv)
   if length(nested)>0
      # nested variable
      for i=1:n
         loc = Int(nested[i])
         if loc < 1 || loc > lv
            error("invalid category code: $(nested[i])")
         end
         Z[i,loc] = x[i]
      end
   else
      # non-nested case
      for i=1:n
         loc = Int(x[i])
         if loc < 1 || loc > lv
            error("invalid category code: $(x[i])")
         end
         Z[i,loc] = 1.0
      end
   end
   return Z
end

function get_design_matrix(X::Matrix{<:Real}; maxlev::Vector{<:Integer}=Int[], cov::Vector{Bool}=Bool[])
   # check
   (n,neff)=size(X)
   if length(maxlev)>0 && length(maxlev)!=neff
      error("argument maxlev size mismatch")
   end
   if length(cov)>0 && length(cov)!=neff
      error("argument cov size mismatch")
   end
   if length(cov)>0
      covar=cov
   else
      covar=zeros(Bool,neff)
   end
   # levels
   if length(maxlev)>0
      if minimum(maxlev)<0
         error("argument maxlev with 0 or negative value")
      end
      lv = maxlev
   else
      lv = zeros(Int,neff)
      for i=1:neff
         if covar[i]
            lv[i] = 1
         else
            lv[i] = maximum(X[:,i])
         end
      end
   end
   # build Y
   ncol = sum(lv)
   Z = zeros(Float64,n,ncol)
   collst = 0
   for i=1:neff
      colfst = collst + 1
      collst = colfst + lv[i] - 1
      Z[:,colfst:collst] .= get_design_matrix(X[:,i],maxlev=lv[i],cov=covar[i])
   end
   return Z
end

"""
     M = directsum(X...)

Returns the direct sum of supplied matrices.
"""
function directsum(X::Matrix...)
    num = length(X)
    mX = Int[ size(x, 1) for x in X ]
    nX = Int[ size(x, 2) for x in X ]
    m = sum(mX)
    n = sum(nX)
    Tv = promote_type(map(x->eltype(x), X)...)

    M =	  zeros(Tv,m,n)
    prev_i = 0
    prev_j = 0
    for i=1:num
       ifst = prev_i + 1
       ilst = ifst + mX[i] - 1
       jfst = prev_j + 1
       jlst = jfst + nX[i] - 1
       M[ifst:ilst,jfst:jlst] = X[i]
       prev_i =		      ilst
       prev_j =		      jlst
    end
    return M
end
