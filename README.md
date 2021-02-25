# AnimalBreedingTools

[![Build Status](https://travis-ci.org/masuday/AnimalBreedingTools.jl.svg?branch=master)](https://travis-ci.org/masuday/AnimalBreedingTools.jl)

## Quick start

This package provides small functions useful for everyday work in animal breeding (or quantitative genetics).

## Example

```
# traditional bending by Hayes and Hill (1981)
B = [
   50.0 48.0
   48.0 40.0
]
bentB = bending(B)

# bending generalized by Jorjani et al. (2003)
V = [
   50.0 48.0
   48.0 40.0
]
Weight = 1 ./ [
   100.0 80.0
    80.0 80.0
]
bentV = bending2(V,Weight)

# normalized Legendre Polynomials
M = [normalized_legendre(x,n) for x=(-1:0.1:1),n=0:3]

# `vech` function
V = [
   50.0 48.0
   48.0 40.0
]
v = vech(V)

# covariance matrix to correlation matrix
# same functionality as StatsBase.cov2cor
R = v2r(bentV)

# correlation matrix to covariance matrix
# same functionality as StatsBase.cor2cov
C = r2v(R,sqrt.(diag(V)))

# direct sum
sV = directsum(V,V)
tV = directsum(V,V,V)
```
