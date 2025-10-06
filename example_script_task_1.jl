# example script task 1
using Pkg
Pkg.instantiate()
Pkg.activate(".")
include("poly_factorization_project.jl")

# DensePolynomials are constructed with lists of Terms, wherein a Term is constructed with a tuple of integers, corresponding to coefficient, degree
f = PolynomialDense([Term(2, 1), Term(1, 0)])
g = PolynomialDense([Term(5, 0), Term(8, 1), Term(3, 2)])
h = PolynomialDense([Term(1, 2)])
println("f(x) = ", f)
println("g(x) = ", g)
println("h(x) = ", h)

# # We can perform the typical operations on these polynomials
println("\nBasic polynomial operations")
println("f(x) plus g(x) = ", f + g)
println("f(x) multiplied by h(x) = ", f * h)
println("g(x) multiplied by h(x) = ", g * h)

# We can take the derivative of these polynomials as well as the derivatives of their products
println("\nDerivative showcase")
f_times_g = f * g
println("The product of f(x) and g(x) = ", f * g)
derivative_f_times_g = derivative(f_times_g)
println("The derivative of f(x)*g(x) = ", derivative_f_times_g)
derivative_prod = (f * derivative(g)) + (derivative(f) * g)
println("The derivative explicitly using the product rule = ", derivative_prod)
println("We can see that ", derivative_f_times_g, " = ", derivative_prod)
println("\n")

#####################
# computing modulus #
#####################

nums = [5, 17, 101]
for num in nums
    f_mod = mod(f, num)
    q, r = div_rem_mod_p(f * h, h, num)
    println("Quotient = ", q)
    println("Remainder = ", r)
    println("f(x) mod $(num) = ", f_mod)
    println("We can see that ", q, " = ", f_mod)
end

println("\nGCD Showcase")
# #############
# # GCD
# #############
nums = [5, 11, 13]
for num in nums
    println("The gcd of ", f * h, " and ", g * h, " mod ", num, " is ", gcd_mod_p(f*h, g*h, num))
end