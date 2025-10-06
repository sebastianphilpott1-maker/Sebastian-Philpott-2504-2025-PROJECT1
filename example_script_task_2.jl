using Pkg
Pkg.instantiate()
Pkg.activate(".")
include("poly_factorization_project.jl")

println("--- BEGIN EXAMPLE SCRIPT 2 ---")

# Example 1: Overflow into negatives by addition of two positive polynomials 
# (2^63 - 1 is the coefficient of f, the largest positive int64 value possible before overflow)
f = PolynomialDense{Int, Int}([Term(9223372036854775807, 2)])
g = PolynomialDense{Int, Int}([Term(1, 2)])
println("We expect that ", f, " + ", g, " is 9223372036854775808x^2")
println("Instead, f + g = ", f + g)
println("\n")

# Example 2: Overflow into positives by addition of two negative polynomials
h = PolynomialDense{Int, Int}([Term(-1, 2)])
i = PolynomialDense{Int, Int}([Term(-9223372036854775808, 2)])
println("We would expect that ", h, " + ", i, " is -9223372036854775809x^2")
println("Instead, h + i = ", h + i)
println("\n")

# Example 3: Overflow into negatives by addition two polynomials 
# whose coefficients are the midpoint between 0 and 2^63 - 1
j = PolynomialDense{Int, Int}([Term(4611686018427387904, 2)])
k = PolynomialDense{Int, Int}([Term(4611686018427387904, 2)])
println("We would epxect that ", j, " + ", k, " is 9223372036854775808x^2")
println("Instead, j + k = ", j + k)
println("\n")

# Example 4: Negative overflow to wrap around close to zero
l = PolynomialDense{Int, Int}([Term(9223372036854775807, 2)])
m = PolynomialDense{Int, Int}([Term(9223372036854775807, 2)])
println("We would epxect that ", l, " + ", m, " is ", BigInt(2) * BigInt(9223372036854775807), "x^2")
println("Instead, l + m = ", l + m)
println("\n")

# Example 5: Positive overflow once more
n = PolynomialDense{Int, Int}([Term(-2305843009213693952, 2)])
o = PolynomialDense{Int, Int}([Term(-9223372036854775807, 2)])
println("We would epxect that ", n, " + ", o, " is ", BigInt(-2305843009213693952) + BigInt(-9223372036854775807), "x^2")
println("Instead, n + o = ", n + o)
println("\n")

println("\nShowcase of correct overflow of polynomial dense using BigInt")
println("All examples use the same numbers as above.")
println("The only difference is they are cast to BigInt\n")

# Examples working with BigInt
f = PolynomialDense{BigInt, Int}([Term(BigInt(9223372036854775807), 2)])
g = PolynomialDense{BigInt, Int}([Term(BigInt(1), 2)])
println("Now we see that ", f, " + ", g, " is ", f + g)
println("\n")


h = PolynomialDense{BigInt, Int}([Term(BigInt(-1), 2)])
i = PolynomialDense{BigInt, Int}([Term(BigInt(-9223372036854775808), 2)])
println("Now we see that ", h, " + ",i, " is ", h + i)
println("\n")

j = PolynomialDense{BigInt, Int}([Term(BigInt(4611686018427387904), 2)])
k = PolynomialDense{BigInt, Int}([Term(BigInt(4611686018427387904), 2)])
println("Now we see that ", j, " + ", k, " is ", j + k)
println("\n")


l = PolynomialDense{BigInt, Int}([Term(BigInt(9223372036854775807), 2)])
m = PolynomialDense{BigInt, Int}([Term(BigInt(9223372036854775807), 2)])
println("Now we see that ", l, " + ", m, " is ", l + m)
println("\n")


n = PolynomialDense{BigInt, Int}([Term(BigInt(-2305843009213693952), 2)])
o = PolynomialDense{BigInt, Int}([Term(BigInt(-9223372036854775807), 2)])
println("Now we see that ", n, " + ", o, " is ", n + o)
println("\n")

println("--- END EXAMPLE SCRIPT 2 ---")