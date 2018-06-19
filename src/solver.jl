"""
   x = solve_pcg(A,b; maxiter=5000,crit=1e-12)

Solve the symmetric positive definite system of equations:
`A` as sparse left-hand side and `b` as a vector of right-hand side.
The Jacobi preconditioner (`diag(A)`) is used.
Convergence criterion `crit` is `norm(b-A*x)/norm(b)`.


```juliadoctest
julia> x = solve_pcg(A,b)
```
"""
function solve_pcg(A, b::Vector{Tv};
                   maxiter=5000,crit=1e-15,verbose=false) where {Tv}
   n = size(A)[1]
   x = zeros(Tv,n,1)
   r = zeros(Tv,n,1)
   p = zeros(Tv,n,1)
   z = zeros(Tv,n,1)
   w = zeros(Tv,n,1)
   Mi = inv(diagm(diag(A)))
   r .= b
   oldtau = 1
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
