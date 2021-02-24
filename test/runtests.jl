using AnimalBreedingTools
using LinearAlgebra
using Printf
using Test

@testset "get_design_matrix" begin
    x = [1,2,1,1,2]
    y = [0.5,0.6,0.7,0.8,0.9]
    # class variable
    X = get_design_matrix(x)
    @test X == [1.0 0.0; 0.0 1.0; 1.0 0.0; 1.0 0.0; 0.0 1.0]
    # covariate
    X = get_design_matrix(x, cov=true)
    @test X == Matrix(convert.(Float64, x[:,:]))
    # nested covariate
    X = get_design_matrix(y, x, cov=true)
    @test X == [0.5 0.0; 0.0 0.6; 0.7 0.0; 0.8 0.0; 0.0 0.9]
end

@testset "Jacobi CG" begin
   n = 10
   niter = 5
   for iter=1:niter
      A = Symmetric(rand(n,n).>0.4) + n*I
      x = ones(n)
      b = A*x
      s = solve_pcg(A, b, eps=1e-16)
      @test x ≈ s
   end
end

@testset "Approximated trace of the inverse" begin
   Ai=[
      7.11933   -2.11262   -1.37931   -1.37931   -0.746269   0.0        0.0        0.0        0.0        0.0        0.0        0.0        0.0      0.0
      -2.11262    6.47343    0.0        0.0       -0.746269  -1.62602   -1.62602    0.0        0.0        0.0        0.0        0.0        0.0      0.0
      -1.37931    0.0        2.06897    0.689655   0.0        0.0        0.0       -1.37931    0.0        0.0        0.0        0.0        0.0      0.0
      -1.37931    0.0        0.689655   2.74237    0.673401   0.0        0.0       -1.37931   -1.3468     0.0        0.0        0.0        0.0      0.0
      -0.746269  -0.746269   0.0        0.673401   3.56661    0.0        0.0        0.673401  -1.3468     0.727273   0.0       -1.45455   -1.3468   0.0
       0.0       -1.62602    0.0        0.0        0.0        2.43902    0.813008   0.0        0.0       -1.62602    0.0        0.0        0.0      0.0
       0.0       -1.62602    0.0        0.0        0.0        0.813008   2.43902    0.0        0.0       -1.62602    0.0        0.0        0.0      0.0
       0.0        0.0       -1.37931   -1.37931    0.673401   0.0        0.0        4.11345    0.681431   0.0       -1.36286    0.0       -1.3468   0.0
       0.0        0.0        0.0       -1.3468    -1.3468     0.0        0.0        0.681431   3.37503    0.0       -1.36286    0.0        0.0      0.0
       0.0        0.0        0.0        0.0        0.727273  -1.62602   -1.62602    0.0        0.0        3.97931    0.0       -1.45455    0.0      0.0
       0.0        0.0        0.0        0.0        0.0        0.0        0.0       -1.36286   -1.36286    0.0        3.47725    0.751527   0.0     -1.50305
       0.0        0.0        0.0        0.0       -1.45455    0.0        0.0        0.0        0.0       -1.45455    0.751527   3.66062    0.0     -1.50305
       0.0        0.0        0.0        0.0       -1.3468     0.0        0.0       -1.3468     0.0        0.0        0.0        0.0        2.6936   0.0
       0.0        0.0        0.0        0.0        0.0        0.0        0.0        0.0        0.0        0.0       -1.50305   -1.50305    0.0      3.00611
   ]
   A = inv(Ai)
   true_tr = sum(diag(A))
   maxiter = 100
   tr = 0
   for i=1:maxiter
      tr = tr + approximate_trace_of_inverse(Ai)
   end
   tr = tr/maxiter
   println("true_tr=$(true_tr)  tr=$(tr)")
   @test abs(true_tr - tr)< 0.1
end

@testset "Symmetric tridiagonal tools" begin
   T = SymTridiagonal([2.0,3.0,6.0,7.0,9.0],[1.0,2.0,2.0,3.0])
   n = size(T,1)
   dv,ev = chol_symtrid(T)
   @test Matrix( cholesky(Matrix(T)).L ) ≈ tril(Matrix(SymTridiagonal(dv,ev)))
   dv,ev = ldlt_symtrid(T)
   D = diagm(0=>dv)
   L = tril(Matrix(SymTridiagonal(ones(n),ev)))
   @test L*D*L' ≈ Matrix(T)
   takahashi_ldlt_symtrid!(dv,ev)
   @test Matrix(SymTridiagonal(Symmetric(inv(Matrix(T))))) ≈ Matrix(SymTridiagonal(dv,ev))
end

@testset "normalized_legendre" begin
   ref=[
      0.70711  -1.22474    1.58114    -1.87083
      0.70711  -1.10227    1.13052    -0.883968
      0.70711  -0.979792   0.727324   -0.149668
      0.70711  -0.857318   0.371568    0.360133
      0.70711  -0.734844   0.0632456   0.673497
      0.70711  -0.61237   -0.197642    0.818486
      0.70711  -0.489896  -0.411096    0.823164
      0.70711  -0.367422  -0.577116    0.715591
      0.70711  -0.244948  -0.695702    0.523831
      0.70711  -0.122474  -0.766853    0.275947
      0.70711   0.0       -0.79057     0.0
      0.70711   0.122474  -0.766853   -0.275947
      0.70711   0.244948  -0.695702   -0.523831
      0.70711   0.367422  -0.577116   -0.715591
      0.70711   0.489896  -0.411096   -0.823164
      0.70711   0.61237   -0.197642   -0.818486
      0.70711   0.734844   0.0632456  -0.673497
      0.70711   0.857318   0.371568   -0.360133
      0.70711   0.979792   0.727324    0.149668
      0.70711   1.10227    1.13052     0.883968
      0.70711   1.22474    1.58114     1.87083
   ]
   leg = [normalized_legendre(x,n) for x=(-1:0.1:1),n=0:3]
   @test isapprox(ref,leg,atol=1e-3)
end
