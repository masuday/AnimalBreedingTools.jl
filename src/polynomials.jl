function normalized_legendre(x,n)
   a = sqrt((2*n+1)/2)
   return a*legendre(x,n)
end
