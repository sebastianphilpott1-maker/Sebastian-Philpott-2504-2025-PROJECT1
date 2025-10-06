#############################################################################
#############################################################################
#
# This file implements operations for polynomials over a field, e.g. Zp[x].
#                                                                               
#############################################################################
#############################################################################

#= TASK LIST:
NOTE - the following functions will NOT work when C is an integer until you 
finish implementing all tasks (and even then, some functions such as gcd simply
cannot be implemented over the ring of integers - see bonus task 3 for details).

TODO (Task 2) 
    In this task, simply edit this file such that the first argument (::Type{C})
    to each function is also the type of the coefficients of the polynomials.
    I.e., we want the type signatures to contain:
        {C, D, P <: Polynomial{C, D}}

TODO (Task 6)
    In this task, override the unimplemented functions (after duplicating this 
    file as per the instructions) for C <: ZModP. When overriding, use the
    corresponding _mod_p function as a guide.

TODO (Task 7)
    Here, (in the duplicated file as per the instructions) override the factor
    function for C <: Integer and implement the Chinese Remainder Theorem (CRT).

=#


"""
Returns the factors of `f` in an array of tuples. 

For a given tuple (g, n) in the returned array, g is an irreducible factor
of f and the multiplicity of g in f is n.

The returned array may also contain a constant for reconstruction of the
leading coefficient.

NOTE: Override this in Task 6 for Zp[x] and in Task 7 for Z[x].
"""
function factor(::Type{C}, f::P)::Vector{Tuple{P, Integer}} where {C, D, P <: Polynomial{C, D}}
    not_implemented_error(f, "factor")
end

"""
Returns the quotient and remainder of num divided by den (where num/den have the 
same concrete type). 

NOTE: Override this in Task 6 for Zp[x].
"""
function div_rem(::Type{C}, num::P, den::P)::Tuple{P, P} where {C, D, P <: Polynomial{C, D}}
    not_implemented_error(num, "div_rem")
end

"""
Distinct degree factorization.

Given a square free polynomial `f` returns a list, `g` such that `g[k]` is a product 
of irreducible polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2, such that the 
product of the list is equal to `f`.

NOTE: Override this in Task 6 for Zp[x].
"""
function dd_factor(::Type{C}, f::P)::Array{P} where {C, D, P <: Polynomial{C, D}}
    not_implemented_error(f, "dd_factor")
end

"""
Distinct degree split.

Returns a list of irreducible polynomials of degree `d` so that the product of
that list is the polynomial `f`.

NOTE: Override this in Task 6 for Zp[x].
"""
function dd_split(::Type{C}, f::P, d::Integer)::Vector{P} where {C, D, P <: Polynomial{C, D}}
    not_implemented_error(f, "dd_split")
end

""" 
Returns the quotient of num divided by den (where num/den have the 
same concrete type) 
"""
div(::Type{C}, num::P, den::P) where {C, D, P <: Polynomial{C, D}} = first(div_rem(C, num, den))

""" 
Returns the remainder of num divided by den (where num/den have the 
same concrete type) 
"""
rem(::Type{C}, num::P, den::P) where {C, D, P <: Polynomial{C, D}} = last(div_rem(C, num, den))

"""
The extended euclid algorithm for polynomials (of the same concrete subtype).

https://edstem.org/au/courses/24848/discussion/2979071

learnt this too late idk if this would fix anything but w/e
"""
function extended_euclid_alg(::Type{C}, f::P, g::P) where {C, D, P <: Polynomial{C, D}}
    return ext_euclid_alg(f, g, rem, div)
end

"""
The greatest common divisor of two polynomials (of the same concrete subtype).
"""
gcd(::Type{C}, f::P, g::P) where {C, D, P <: Polynomial{C, D}} = extended_euclid_alg(C, f, g) |> first

"""
Yun's algorithm to compute a square free polynomial can be performed over any so-called
`perfect field`. All the usual fields you would have worked with are perfect. In particular,
any field with characteristic zero, or finite fields (i.e., Zp) are perfect fields. Thus,
we apply Yun's algorithm.

Precondition: C represents a perfect field.

Note: to account for the case of fields with characteristic p (e.g., Zp), we first divide
by the minimum degree to remove extraneous factors of x. These can cause issues due to the 
following lemma:

Lemma:
Given F a characteristic zero or finite field and f a polynomial in F[x],
    1) If char(F) = 0 then: f' = 0 implies f is constant
    2) If char(F) > 0 then: f' = 0 implies f has a factor of x^char(F)

https://en.wikipedia.org/wiki/Square-free_polynomial
https://en.wikipedia.org/wiki/Perfect_field
"""
function square_free(::Type{C}, f::P) where {C, D, P <: Polynomial{C, D}}
    # Remove minimum degree (in case char(C) != 0)
    min_deg = last(f).degree
    vt = filter(t -> !iszero(t), collect(f))
    f = P( map(t -> Term(t.coeff, t.degree - min_deg), vt) )

    # Compute the gcd of f, f'
    der_f = derivative(f)
    sqr_part = gcd(C, f, der_f)

    iszero(sqr_part) && return f * (min_deg > zero(min_deg) ? x_poly(P) : one(P))

    # Remove factors with multiplicity > 1
    sqr_free = div(f, sqr_part)

    # Add one factor of `x` back in if necessary
    if min_deg > zero(min_deg) 
        sqr_free *= x_poly(P)
    end

    return sqr_free
end

"""
Compute the number of times g divides f.
"""
function multiplicity(::Type{C}, f::P, g::P)::Integer where {C, D, P <: Polynomial{C, D}}
    degree(gcd(C, f, g)) == 0 && return 0
    return 1 + multiplicity(C, div(f, g), g)
end
