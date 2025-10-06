#############################################################################
#############################################################################
#
# This file implements polynomial division for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""  Modular algorithm.
f divide by g

f = q*g + r

p is a prime

Returns quotient and remainder (q, r)

Precondition: 
    1) degree(num) == degree(mod(num, prime))
    2) mod(den, prime) != 0

Note: clearly a monic polynomial will satisfy 1) in Zp for any field.
"""
function div_rem_mod_p(num::P, den::P, prime::Integer)::Tuple{P, P} where {C, D, P <: Polynomial{C, D}}
    if C <: ZModP
        return div_rem(C, num, den)
    else
        num_poly = convert_polynomial(Val(prime), num)
        den_poly = convert_polynomial(Val(prime), den)
        q, f = div_rem(ZModP{C, prime}, num_poly, den_poly)
        q = convert_polynomial(Val(prime), q)
        f = convert_polynomial(Val(prime), f)
        return q, f
    end
end

"""
The quotient from polynomial division modulo a prime. 
"""
div_mod_p(num::P, den::P, prime::Integer) where {C, D, P <: Polynomial{C, D}} = first(div_rem_mod_p(num, den, prime))

"""
The remainder from polynomial division modulo a prime.
"""
rem_mod_p(num::P, den::P, prime::Integer) where {C, D, P <: Polynomial{C, D}} = last(div_rem_mod_p(num, den, prime))
