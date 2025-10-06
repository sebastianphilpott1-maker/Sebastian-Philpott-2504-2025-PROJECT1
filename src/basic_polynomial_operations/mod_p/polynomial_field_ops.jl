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
This feels like an easy way to lose marks but I'd be remiss if I didn't say this.

I do not really understand the maths in here. I understand modulus.
Beyond that, I'm lost. Most of the code is copied from the mod_p versions.
"""


"""
Returns the factors of `f` in an array of tuples. 

For a given tuple (g, n) in the returned array, g is an irreducible factor
of f and the multiplicity of g in f is n.

The returned array may also contain a constant for reconstruction of the
leading coefficient.

NOTE: Override this in Task 6 for Zp[x] and in Task 7 for Z[x].
"""
function factor(::Type{C}, f::P)::Vector{Tuple{P, Integer}} where {T <: Integer, N, D, C <: ZModP{T, N}, P <: Polynomial{C, D}}
    # Cantor Zassenhaus factorization
    # REMEMBER WE ARE ALWAYS WORKING IN MOD N. ANY OPS WE DO ON ZMODP
    # OBJECTS HAVE MOD OPS BAKED INTO IT. DONKEY DELETE MOD OPS HERE
    degree(f) ≤ 1 && return [(f, 1)]

    ret_val = Tuple{P, Int}[]

    # make f square-free
    f_sf = square_free(C, f)

    # Drop to Zp[x] - WE ARE ALREADY IN ZP CHIEF 

    # make f monic - USE CURRENT REDUCED f_copy NOT f
    old_coeff = leading(f_sf).coeff # ZMODP OBJECT DONT FORGET
    f_sfm = f_sf * inv(old_coeff)

    dds = dd_factor(C, f_sfm)

    for (k, dd) in enumerate(dds)
        sp = dd_split(C, dd, k)
        sp = map((p) -> p*inv(leading(p).coeff), sp)
        for mp in sp
            push!(ret_val, (mp, multiplicity(C, f_sfm, mp)) )
        end
    end
    
    # Append the leading coefficient as well
    push!(ret_val, (leading(f).coeff * one(f), 1) )

    return ret_val
end

"""
Returns the quotient and remainder of num divided by den (where num/den have the 
same concrete type). 

NOTE: Override this in Task 6 for Zp[x].
"""
function div_rem(::Type{C}, num::P, den::P)::Tuple{P, P} where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    f = num
    g = den

    # If we assume num & den are polynomials with ZModP coeffs, no mod needed
    # nor is the assertion since by definition num is already mod N
    iszero(g) && throw(DivideError())
    iszero(f) && return zero(P), zero(P)

    q = zero(P)
    prev_degree = degree(f)
    while !iszero(f) && degree(f) ≥ degree(g)
        # im so disgustingly lost

        # THESE SHOULD BE TERMS
        # @show typeof(leading(f))
        # @show typeof(leading(g))

        # THIS TOO SHOULD BE A TERM, CALLING TERM DIV ON TWO POLYS WITH 
        # ZMODP COEFFICIENTS.
        term = div(leading(f), leading(g))
        
        
        # TY CONVENIENCE CONSTRUCTOR
        h = P([term])  #syzergy
        f = f - h*g
        q = q + h
        
        
        # we are working within mod N already so just do the operations
        # that are defined for zmodp objects and all the modulus stuff
        # is handled
        prev_degree == degree(f) && break
        prev_degree = degree(f)

    end
    @assert iszero(num  - (q*g + f))
    return q, f
end

"""
Distinct degree factorization.

Given a square free polynomial `f` returns a list, `g` such that `g[k]` is a product 
of irreducible polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2, such that the 
product of the list is equal to `f`.

NOTE: Override this in Task 6 for Zp[x].

Anti-dizzy Note: In our g, as per lecture slides, we are collecting the g_i

We can optimise the loop by stopping when i > deg(f) / 2 since by that point
we will be looking at duplicates. So, if we have square free monic poly and we
know that ddf decomposes polys into groups of products of irreducible poly of
equiv degree

https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Distinct-degree_factorization

ALSO also, looking over everything before submitting, pretty sure this is awful. Like all of it
is just horrid and doesn't work and I don't get it. Removing pow_mod was a mistake. Submit first, try fix later
"""
function dd_factor(::Type{C}, f::P)::Array{P} where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    x = x_poly(P)
    w = deepcopy(x)
    g = []
    max_deg = degree(f)
    f_copy = f

    #Looping over degrees
    k = (max_deg) ÷ 2
    for i in 1:k
        # what have i done where is pow_mod
        w = rem(C, w^N, f_copy)
        push!(g, gcd(C, f_copy, w - x))
        f_copy = div(C, f_copy, g[i])
    end

    #edge case for final factor
    f_copy != one(P) && push!(g, f_copy)
    
    return g
end


"""
Distinct degree split.

Returns a list of irreducible polynomials of degree `d` so that the product of
that list is the polynomial `f`.

NOTE: Override this in Task 6 for Zp[x].
"""
function dd_split(::Type{C}, f::P, d::Int)::Vector{P} where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    # absolutely NO idea if this works correctly or not
    degree(f) == D(d) && return [f]
    degree(f) == D(0) && return P[]

    w = rand(P, degree = d, monic = true)
    power = (N^d - 1) ÷ 2

    # because original things get modded and we use pow mod but we dont have
    # it defined for polynomials ^ polynomials, we do it here could use helper
    # but im over this
    dummy_power = power
    a = one(P)
    base = rem(C, w, f)
    while dummy_power > 0
        if (dummy_power & 1) == 1
            a = rem(C, a*base, f)
        end
        dummy_power >>= 1
        base = rem(C, base*base, f)
    end

    # gcd on zmodp objects has mod built in
    g = gcd(C, a - one(P), f)
    
    gbar = div(C, f, g)
    return vcat(dd_split(C, g, d), dd_split(C, gbar, d))
end

""" 
Returns the quotient of num divided by den (where num/den have the 
same concrete type) 
"""
div(::Type{C}, num::P, den::P) where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}} = first(div_rem(C, num, den))

""" 
Returns the remainder of num divided by den (where num/den have the 
same concrete type) 
"""
rem(::Type{C}, num::P, den::P) where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}} = last(div_rem(C, num, den))

"""
The extended euclid algorithm for polynomials (of the same concrete subtype).

these people really seem to care about explaining whats going on pretty concisely
https://discourse.julialang.org/t/can-someone-explain-closures-to-me/105605
"""
function extended_euclid_alg(::Type{C}, f::P, g::P) where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    # im still not sure if i understand closure but rem_function and
    # div function calls within ext_euclid are two parameter functions

    # but our rem and div (not the mod_p versions) are three params
    # ENTER CLOSURE 
    # we do below so that when we map dummy_rem and dummy_div to rem and div
    # dummy versions store the reference the variable C from the outer
    # scope. So, when we call ext_euclid and it calls div and rem_function
    # since the reference to C is stored, it infers which div and rem to use
    # namely, the ones defined in here

    dummy_div = (a, b) -> div(C, a, b)
    dummy_rem = (a, b) -> rem(C, a, b)
    return ext_euclid_alg(f, g, dummy_rem, dummy_div)
end

"""
The greatest common divisor of two polynomials (of the same concrete subtype).
"""
gcd(::Type{C}, f::P, g::P) where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}} = extended_euclid_alg(C, f, g) |> first

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
function square_free(::Type{C}, f::P) where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    # Remove minimum degree (in case char(C) != 0)
    min_deg = last(f).degree
    vt = filter(t -> !iszero(t), collect(f))
    f = P( map(t -> Term(t.coeff, t.degree - min_deg), vt) )

    # Compute the gcd of f, f'
    der_f = derivative(f)
    sqr_part = gcd(C, f, der_f)

    iszero(sqr_part) && return f * (min_deg > zero(min_deg) ? x_poly(P) : one(P))

    # Remove factors with multiplicity > 1
    sqr_free = div(C, f, sqr_part)

    # Add one factor of `x` back in if necessary
    if min_deg > zero(min_deg) 
        sqr_free *= x_poly(P)
    end

    return sqr_free
end

"""
Compute the number of times g divides f.
"""
function multiplicity(::Type{C}, f::P, g::P)::Integer where {T <: Integer, N, C <: ZModP{T, N}, D, P <: Polynomial{C, D}}
    degree(gcd(C, f, g)) == 0 && return 0
    return 1 + multiplicity(C, div(C, f, g), g)
end
