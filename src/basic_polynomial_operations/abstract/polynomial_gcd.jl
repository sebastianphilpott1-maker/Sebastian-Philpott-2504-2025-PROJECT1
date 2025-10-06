#############################################################################
#############################################################################
#
# This file implements polynomial GCD for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
The extended euclid algorithm for polynomials (of the same concrete subtype) modulo prime.

Note, when working with polynomials over Z[x] we cannot compute the extended euclidean algorithm (EEA) directly.
This is because we do not have division in Z[x]. Hence, we must drop to the field Zp[x] and compute the EEA
with respect to some prime. 
"""
function extended_euclid_alg_mod_p(a::P, b::P, prime::Integer) where {C, D, P <: Polynomial{C, D}}
    # if we have a polynomial that is already zmodp obj coeff, don't try convert it obviously
    if C <: ZModP
        return extended_euclid_alg(C, a, b)
    else
        # convert polynomials
        poly_a = convert_polynomial(Val(prime), a)
        poly_b = convert_polynomial(Val(prime), b)

        # call the one for zmodp objects
        g, s, t = extended_euclid_alg(ZModP{C, prime}, poly_a, poly_b)

        # convert back
        g = convert_polynomial(Val(prime), g)
        s = convert_polynomial(Val(prime), s)
        t = convert_polynomial(Val(prime), t)
        return g, s, t
    end
end

"""
The GCD of two polynomials (of the same concrete subtype) modulo prime.
"""
gcd_mod_p(a::P, b::P, prime::Integer) where {C, D, P <: Polynomial{C, D}} = extended_euclid_alg_mod_p(a, b, prime) |> first
