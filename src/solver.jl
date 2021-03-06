"""
    x = solve_pcg(A,b; maxiter=5000,eps=1e-15,reset=50,verbose=false)

Solve the symmetric positive definite system of equations:
`A` as sparse left-hand side and `b` as a vector of right-hand side.
The Jacobi preconditioner (`diag(A)`) is used.
The maximum number of iterations is specified with `maxiter`.
Convergence criterion `eps` is `norm(b-A*x)/norm(b)`.
The residual vector will be recalculated every `reset` rounds.
The details will be shown with `verbose=true`.

```juliadoctest
julia> x = solve_pcg(A,b)
```
"""
function solve_pcg(A, b::Vector{Tv};
                   maxiter=5000,eps=1e-15,reset=50,verbose=false) where {Tv}
   n = size(A,1)
   x = zeros(Tv,n)
   r = zeros(Tv,n)
   p = zeros(Tv,n)
   z = zeros(Tv,n)
   w = zeros(Tv,n)
   Mi = Diagonal(1 ./ diag(A))
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
      if i % reset == 0
         r .= b - A*x
      else
         r .= r - alpha*w
      end
      conv = norm(r)^2/norm(b)^2
      if verbose
         print(@sprintf("round %5d: conv %8.3e\n",i,conv))
      end
      if conv<eps; break; end
      oldtau = tau
   end
   x
end

"""
   tr = approximate_trace_of_inverse(A)
"""
function approximate_trace_of_inverse(A; rounds=1000, msglev=0)
   if rounds<1
      throw( ArgumentError("negative rounds") )
   end
   tr = Threads.Atomic{Float64}(0.0)
   done = Threads.Atomic{Int}(0)
   n = size(A,1)
   @threads for iter=1:rounds
      atomic_add!(done, 1)
      z = [ifelse(rand()>0.5,1.0,-1.0) for i in 1:n]
      #tr = tr + z'*y
      if msglev>0
         @time y = solve_pcg(A,z)
         println("iter $(done.value) thread $(Threads.threadid())")
      else
         y = solve_pcg(A,z)
      end
      atomic_add!(tr, z'*y)
   end
   return tr.value/rounds
end
