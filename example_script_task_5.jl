using Pkg
Pkg.activate(".")
include("poly_factorization_project.jl")
# am i dumb or is there no script 4??????????

println("---BEGIN ZMODP OBJECT OPERATIONS TESTS---")

# Define two zmodp objects
a = ZModP{BigInt, 7}(BigInt(11))
b = ZModP{BigInt, 7}(BigInt(17))

# Addition of two zmodp objects and eq
@assert (a + b) == mod(11 + 17, 7)
@assert (b + a) == mod(17 + 11, 7)
@assert (a + b) == (b + a)
println("11 mod 7 + 17 mod 7 is ", a + b, " which is the same as mod(11 + 17, 7) = ", mod(11 + 17, 7))

# Addition of a zmodp object and integer
@assert (a + 17) == mod(11 + 17, 7)
@assert (11 + b) == mod(11 + 17, 7)

# Subtraction operations
@assert (-a) == mod(-11, 7)
@assert (a - b) == mod(11 - 17, 7)
@assert (b - a) == mod(17 - 11, 7)
@assert (a - 17) == mod(11 - 17, 7)
@assert (11 - b) == mod(11 - 17, 7)
println("11 mod 7 - 17 mod 7 is ", a - b, " which is the same as mod(11 - 17, 7) = ", mod(11 - 17, 7))

# Multiplication between two zmodp objects
@assert (a * b) == mod(11 * 17, 7)
@assert (b * a) == mod(17 * 11, 7)
@assert (a * b) == (b * a)
println("11 mod 7 * 17 mod 7 is ", a * b, " which is the same as mod(11 * 17, 7) = ", mod(11 * 17, 7))

# Multiplication between int and zmodp
@assert (a * 17) == mod(11 * 17, 7)
@assert (11 * b) == mod(11 * 17, 7)

# inverse - Now that I'm testing it after everything is nearly done
# I don't think this works correctly..?
# this is meant to yield a^{-1}
@assert inv(a) == mod(int_inverse_mod(11, 7), 7)

# And this should yield 1 itself
@assert get_zmod_val(a * inv(a)) == mod(11 * int_inverse_mod(11, 7), 7)

# Division between two zmodp objects
@assert (a ÷ b) == mod(11 * int_inverse_mod(17, 7), 7)
@assert (b ÷ a) == mod(17 * int_inverse_mod(11, 7), 7)
println("11 mod 7 ÷ 17 mod 7 is ", a ÷ b, " which is the same as 11 mod 7 * inverse of 17 mod 7 = ", mod(11 * int_inverse_mod(17, 7), 7))

# Division between zmodp obj and integer
@assert (a ÷ 17) == mod(11 * int_inverse_mod(17, 7), 7)
@assert (17 ÷ a) == mod(17 * int_inverse_mod(11, 7), 7)

# Crude test number 1 of ^
@assert a*a*a*a*a == a^5
println("11 mod 7 to the power of 5 is ", a^5, " which is the same as a*a*a*a*a = ", a*a*a*a*a)

# Less crude test but still crude
@assert a^0 == one(a)
@assert a^1 == a
@assert (a^2) * (a^2) == (a^4)
@assert (a^4)^2 == (a^(4*2))
@assert a^(-3) == (a^(-1))*(a^(-1))*(a^(-1))


