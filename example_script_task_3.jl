using Pkg
Pkg.activate(".")
include("poly_factorization_project.jl")

using BenchmarkTools
using Random
using Test

println("TASK 3.7 TO 3.9")

println("---BEGIN INCORRECT SPARSE OVERFLOW DEMONSTRATION---")
println("---BELOW EQUALITIES ASSUME OVERFLOW\n---")

println("--- BEGIN EXAMPLE SCRIPT 2 ---")

# Example 1: Overflow into negatives by addition of two positive polynomials 
# (2^63 - 1 is the coefficient of f, the largest positive int64 value possible before overflow)
f = PolynomialSparse{Int, Int}([Term(9223372036854775807, 2)])
g = PolynomialSparse{Int, Int}([Term(1, 2)])
println("We expect that ", f, " + ", g, " is 9223372036854775808x^2")
println("Instead, f + g = ", f + g)
println("\n")

# Example 2: Overflow into positives by addition of two negative polynomials
h = PolynomialSparse{Int, Int}([Term(-1, 2)])
i = PolynomialSparse{Int, Int}([Term(-9223372036854775808, 2)])
println("We would expect that ", h, " + ", i, " is -9223372036854775809x^2")
println("Instead, h + i = ", h + i)
println("\n")

# Example 3: Overflow into negatives by addition two polynomials 
# whose coefficients are the midpoint between 0 and 2^63 - 1
j = PolynomialSparse{Int, Int}([Term(4611686018427387904, 2)])
k = PolynomialSparse{Int, Int}([Term(4611686018427387904, 2)])
println("We would epxect that ", j, " + ", k, " is 9223372036854775808x^2")
println("Instead, j + k = ", j + k)
println("\n")

# Example 4: Negative overflow to wrap around close to zero
l = PolynomialSparse{Int, Int}([Term(9223372036854775807, 2)])
m = PolynomialSparse{Int, Int}([Term(9223372036854775807, 2)])
println("We would epxect that ", l, " + ", m, " is ", BigInt(2) * BigInt(9223372036854775807), "x^2")
println("Instead, l + m = ", l + m)
println("\n")

# Example 5: Positive overflow once more
n = PolynomialSparse{Int, Int}([Term(-2305843009213693952, 2)])
o = PolynomialSparse{Int, Int}([Term(-9223372036854775807, 2)])
println("We would epxect that ", n, " + ", o, " is ", BigInt(-2305843009213693952) + BigInt(-9223372036854775807), "x^2")
println("Instead, n + o = ", n + o)
println("\n")

println("\nShowcase of correct overflow of polynomial dense using BigInt")
println("All examples use the same numbers as above.")
println("The only difference is they are cast to BigInt\n")

# Examples working with BigInt
f = PolynomialSparse{BigInt, Int}([Term(BigInt(9223372036854775807), 2)])
g = PolynomialSparse{BigInt, Int}([Term(BigInt(1), 2)])
println("Now we see that ", f, " + ", g, " is ", f + g)
println("\n")


h = PolynomialSparse{BigInt, Int}([Term(BigInt(-1), 2)])
i = PolynomialSparse{BigInt, Int}([Term(BigInt(-9223372036854775808), 2)])
println("Now we see that ", h, " + ",i, " is ", h + i)
println("\n")

j = PolynomialSparse{BigInt, Int}([Term(BigInt(4611686018427387904), 2)])
k = PolynomialSparse{BigInt, Int}([Term(BigInt(4611686018427387904), 2)])
println("Now we see that ", j, " + ", k, " is ", j + k)
println("\n")


l = PolynomialSparse{BigInt, Int}([Term(BigInt(9223372036854775807), 2)])
m = PolynomialSparse{BigInt, Int}([Term(BigInt(9223372036854775807), 2)])
println("Now we see that ", l, " + ", m, " is ", l + m)
println("\n")


n = PolynomialSparse{BigInt, Int}([Term(BigInt(-2305843009213693952), 2)])
o = PolynomialSparse{BigInt, Int}([Term(BigInt(-9223372036854775807), 2)])
println("Now we see that ", n, " + ", o, " is ", n + o)
println("\n")

println("--- END EXAMPLE SCRIPT 2 ---")

println("---END CORRECT SPARSE OVERFLOW DEMONSTRATION\n---")

println("---BEGIN BENCHMARK OF SPARSE AND DENSE POLYNOMIALS")
Random.seed!(2504)
max_degree = 100
num_terms = 90

# Build dense and sparse polynomials - make em kinda big but kinda not
polysparse_a = rand(PolynomialSparse{BigInt,Int}; degree = max_degree, terms = num_terms, max_coeff = 100)
polysparse_b = rand(PolynomialSparse{BigInt,Int}; degree = max_degree, terms = num_terms, max_coeff = 100)
polydense_a = rand(PolynomialDense{BigInt,Int}; degree = max_degree, terms = num_terms, max_coeff = 100)
polydense_b = rand(PolynomialDense{BigInt,Int}; degree = max_degree, terms = num_terms, max_coeff = 100)

println("Addition for sparse polynomials of degree $(max_degree) with $(num_terms) terms:")
@btime polysparse_a + polysparse_b

println("Addition for dense polynomials of degree $(max_degree) with $(num_terms) terms:")
@btime polydense_a+ polydense_b

println("Multiplication for sparse polynomials of degree $(max_degree) with $(num_terms) terms:")
@btime polysparse_a * polysparse_b

println("Multiplication for dense polynomials of degree $(max_degree) with $(num_terms) terms:")
@btime polydense_a * polydense_b
println()