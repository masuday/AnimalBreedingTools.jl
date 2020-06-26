"""
   x = solve_pcg(A,b; maxiter=5000,crit=1e-12,verbose=false)

Solve the symmetric positive definite system of equations:
`A` as sparse left-hand side and `b` as a vector of right-hand side.
The Jacobi preconditioner (`diag(A)`) is used.
The maximum number of iterations is specified with `maxiter`.
Convergence criterion `crit` is `norm(b-A*x)/norm(b)`.
The details will be shown with `verbose=true`.


```juliadoctest
julia> x = solve_pcg(A,b)
```
"""
function solve_pcg(A, b::Vector{Tv};
                   maxiter::Int=5000,crit::Float64=1e-15,verbose::Bool=false) where {Tv}
   n = size(A)[1]
   x = zeros(Tv,n)
   r = zeros(Tv,n)
   p = zeros(Tv,n)
   z = zeros(Tv,n)
   w = zeros(Tv,n)
   Mi = Diagonal(1 ./ diag(A))
   size(A),Tv
   r .= b
   oldtau = 1.0
   for i=1:maxiter
      z .= Mi*r
      tau = z'*r
      if i==1
         beta = 0
         p .= z
      else
         beta = tau/oldtau
         p .= z + beta*p
      end
      w .= A*p
      alpha = tau/(p'*w)
      x .= x + alpha*p
      if i%50 == 0
         r .= b - A*x
      else
         r .= r - alpha*w
      end
      conv = norm(r)^2/norm(b)^2
      if verbose
         print(@sprintf("round %5d: conv %8.3e\n",i,conv))
      end
      if conv<crit; break; end
      oldtau = tau
   end
   x
end
