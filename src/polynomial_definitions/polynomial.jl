##################################################################################
##################################################################################
#
# This file defines the abstract polynomial type alongside a number of operations 
# operations that concrete subtypes should implement.
#                                                                               
##################################################################################
##################################################################################

###################
# Polynomial type #
###################

"""
An abstract `Polynomial` type, acting as the supertype for all implementations of a polynomial.


Introduction:
The functions defined in this file, and in the folder src/basic_polynomial_operations/abstract
provide an interface for any subtype of Polynomial. Some operations (e.g. *) can 
be implemented at the abstract level assuming other functions (e.g. iterate) are implemented for
the concrete subtypes of Polynomial. This utilises Julia's ability to do multiple dispatch, where
the compiler can determine the correct function/method to call depending on the type of polynomial
it is working with.

For any function that is not implemented at the abstract level, see the PolynomialDense struct
for a particular implementation.

Note: although a function may be implemented correctly at the abstract level, without leveraging 
the details of the specific implementation the function may be very slow. Consider overriding these
functions for your concrete subtypes.


Constructors:
We assume that at least the following two constructors exist for a given concrete subtype:
    `Polynomial()` -> the empty constructor should produce the zero polynomial
    `Polynomial(vt::Vector{Term})` -> given a vector of terms `vt`, constructs a polynomial with those terms

We assume also that for any concrete subtype of Polynomial, no zero term is stored with degree higher
than the leading term (ie, we do not store 1 + 2x + x^2 + 0x^3 + 0x^7)
"""
abstract type Polynomial{C, D} end

"""
This function maintains the invariant of the Polynomial type so that there are no zero terms beyond the highest
non-zero term.
"""
# my own note: C and D are type variables, and this function works
# for any concrete values of them
# julia will then sub them with Int, UIntx etc whenever called
# Tells the compiler what C and D are supposed to be
function trim!(p::Polynomial{C, D})::Polynomial{C, D} where {C, D}
    # We stop if the leading term has a non-zero coefficient or the polynomial is the zero polynomial
    while iszero(leading(p)) && !iszero(p)
        pop!(p)
    end
    return p
end

#########################
# Abstract Constructors #
#########################

"""
Construct a polynomial with a single term.
"""
function (::Type{P})(t::Term{C, D})::P where {C, D, P <: Polynomial{C, D}}
    return P([t])
end

"""
Construct a polynomial of the form x^p-x.
"""
function cyclotonic_polynomial(::Type{P}, p::D)::P where {C, D, P <: Polynomial{C, D}}
    return P([Term{C, D}(one(C), p), Term{C, D}(-one(C), zero(D))])
end

"""
Construct a polynomial of the form x-n.
"""
function linear_monic_polynomial(::Type{P}, n::C)::P where {C, D, P <: Polynomial{C, D}} 
    return P([Term{C, D}(one(C), one(D)), Term{C, D}(-C(n), zero(D))])
end

"""
Construct a polynomial of the form x.
"""
function x_poly(::Type{P})::P where {C, D, P <: Polynomial{C, D}} 
    return P([Term{C, D}(one(C), one(D))])
end
x_poly(p::P) where {C, D, P <: Polynomial{C, D}} = x_poly(P)

"""
Creates the zero polynomial.
"""
function zero(::Type{P})::P where {C, D, P <: Polynomial{C, D}} 
    return P()
end
zero(p::P) where {C, D, P <: Polynomial{C, D}} = zero(P)

"""
Creates the unit polynomial.
"""
function one(::Type{P})::P where {C, D, P <: Polynomial{C, D}} 
    return P([Term{C, D}(one(C), zero(D))])
end
one(p::P) where {C, D, P <: Polynomial{C, D}} = one(P)

"""
Generates a random polynomial.
"""
function rand(::Type{P} ; 
                degree::Int = -1, 
                terms::Int = -1, 
                max_coeff::Int = 100, 
                mean_degree::Float64 = 5.0,
                prob_term::Float64  = 0.7,
                monic = false,
                condition = (p)->true) where {C, D, P <: Polynomial{C, D}}
        
    while true 
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = rand(1:max_coeff,_terms+1)
        monic && (coeffs[end] = 1)
        p = P( [Term{C, D}(C(coeffs[i]), D(degrees[i])) for i in 1:length(degrees)] )
        condition(p) && return p
    end
end

###########
# Display #
###########

"""
Show a polynomial.
"""
function show(io::IO, p::Polynomial{C, D}) where {C, D}
    # if it is a zero polynomial, show zero
    if iszero(p)
        print(io, "0")
        return
    end

    # sorting based off of degree
    terms_list = [term for term in p]
    sorted_terms = sort(terms_list; by = t -> t.degree, rev = true)

    # flag to check for first term so we know to print - or not
    is_first = true
    for term in sorted_terms
        # https://stackoverflow.com/questions/43577892/finding-out-types-of-variables-in-julia
        # If we have a zmodp coefficient, get its integer
        if isa(term.coeff, ZModP)
            temp_coeff = get_zmod_val(term.coeff)
        else
            temp_coeff = term.coeff
        end

        # check zero coefficients
        if temp_coeff == 0
            continue
        end
        
        # Handle sign printing
        # if it is the first term print -
        if is_first
            if temp_coeff < 0
                print(io, "-")
            end
            is_first = false
        
        # if not first, decide on signage
        else
            if temp_coeff > 0
                print(io, " + ")
            else
                print(io, " - ")
            end
        end

        # take abs since we have already handled - sign
        to_print = abs(temp_coeff)

        # if degree is zero, print constant
        if term.degree == D(0)
            print(io, 1)
        
        # if degree is one, print term
        elseif term.degree == D(1)
            if to_print == 1
                print(io, "x")
            else
                print(io, "$(to_print)x")
            end
        
        elseif to_print == 1
            print(io, "x^$(term.degree)")
        # print all other terms
        else
            print(io, "$(to_print)x^$(term.degree)")
        end
    end         
end

######################
# CONVERT POLYNOMIAL #
######################
"""
Function as described in 6.4

Passed a polynomial of any type, reads the coefficient.

If the coefficient, C, is of type Int or BigInt, returns a polynomial
of the form P{ZModP{C, N}, D}

If the coefficient, C, is of type ZModP{C, N}, returns a polynomial of
the form P{C, D}

I could not implement this a better way I am sorry.
"""
function convert_polynomial(::Val{N}, p::P) where {C, D, P <: Polynomial{C, D}, N}
    # If the polynomial is dense
    if isa(p, PolynomialDense)

        # if the dense poly has terms of ZModP
        if isa(p.terms[1].coeff, ZModP)
            T, blah = get_zmodp_types(p.terms[1].coeff)

            # build a list of terms
            new_terms = Term{T, D}[]
            for term in p.terms
                push!(new_terms, Term{T, D}(get_zmod_val(term.coeff), term.degree))
            end

            # pass this list to the constructor with new params
            return PolynomialDense{T, D}(new_terms)
        
        # if the dense poly has integer coeffs
        else
            new_terms = Term{ZModP{C, N}, D}[]
            for term in p.terms
                push!(new_terms, Term{ZModP{C, N}, D}(ZModP{C, N}(C(term.coeff)), term.degree))
            end
            return PolynomialDense{ZModP{C, N}, D}(new_terms)
        end
    
    # If the polynomial is sparse
    elseif isa(p, PolynomialSparse)

        # if the sparse poly has terms with zmodp coeffs
        if isa(leading(p).coeff, ZModP)
            T, blahblah = get_zmodp_types(leading(p).coeff)
            new_terms = Term{T, D}[]
            for t in p
                push!(new_terms, Term{T, D}(get_zmod_val(t.coeff), t.degree))
            end
            return PolynomialSparse{T, D}(new_terms)
        
        # if it has integer coeffs
        else
            new_terms = Term{ZModP{C, N}, D}[]
            for t in p
                push!(new_terms, Term{ZModP{C, N}, D}(ZModP{C, N}(C(t.coeff)), t.degree))
            end
            return PolynomialSparse{ZModP{C, N}, D}(new_terms)
        end
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the (non-zero) terms of the polynomial in ascending order. 
This implements the iteration interface.

This must be overridden by concrete subtypes.
"""
function iterate(p::Polynomial{C, D}, state=1) where {C, D}
    not_implemented_error(p, "iterate")
end

##############################
# Queries about a polynomial #
##############################

"""
The number of (non-zero) terms of the polynomial.

This must be overridden by concrete subtypes.
"""
length(p::Polynomial{C, D}) where {C, D} = not_implemented_error(p, "length")

"""
The term of smallest degree in this polynomial.

This must be overridden by concrete subtypes.
"""
last(p::Polynomial{C, D}) where {C, D} = not_implemented_error(p, "last")

"""
The leading term of the polynomial.

This must be overridden by concrete subtypes.
"""
leading(p::Polynomial{C, D}) where {C, D} = not_implemented_error(p, "leading")

"""
Returns the coefficients of the polynomial.
"""
function coeffs(p::Polynomial{C, D})::Vector{C} where {C, D}
    [t.coeff for t in p]
end

"""
The degree of the polynomial.
"""
function degree(p::Polynomial{C, D})::D where {C, D}
    leading(p).degree 
end

"""
The content of the polynomial is the GCD of its coefficients.
"""
function content(p::Polynomial{C, D})::C where {C, D} 
    euclid_alg(coeffs(p))
end

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::Polynomial{C, D}, x::T) where {C, D, T <: Number} = sum(evaluate(t, x) for t in f)

################################
# Pushing and popping of terms #
################################

"""
Push a new leading term into the polynomial (note - a constant can be pushed onto the zero polynomial).
This should modify the existing polynomial p in place (ie, no new polynomials should be created).

Note: you may wish to throw an error if pushing a term of degree that is already in the 
polynomial. However, the exact implementation details are up to you.

This must be overridden by concrete subtypes.
"""
function push!(p::Polynomial{C, D}, t::Term{C, D}) where {C, D}
    not_implemented_error(p, "push!")
end

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.

Note - this should modify the existing polynomial p in place (ie, no new polynomials should be 
created).

This must be overridden by concrete subtypes.
"""
function pop!(p::Polynomial{C,D})::Term{C, D} where {C, D} 
    not_implemented_error(p, "pop")
end

"""
Check if the polynomial is zero.
"""
function iszero(p::Polynomial{C, D})::Bool where {C, D}
    iszero(leading(p)) && (degree(p) == D(0))
end

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
function -(p::P) where {C, D, P <: Polynomial{C, D}} 
    P(map((pt)->-pt, p))
end

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::P)::P where {C, D, P <: Polynomial{C, D}} 
    der_p = P()
    for term in p
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p, der_term)
    end
    return trim!(der_p)
end

"""
Returns a primitive polynomial.
E.g., 
    f = 6x^3 + 6x^2 + 9x + 12
    prim_part(f) = 2x^3 + 2x^2 + 3x + 4
"""
function prim_part(f::P) where {C, D, P <: Polynomial{C, D}}
    iszero(f) && return f
    content = gcd(map(t -> t.coeff, f))
    return P( map(t -> Term{C, D}(C(t.coeff ÷ content), D(t.degree)), f) ) # Exact integer division
end

"""
A square free polynomial modulo a prime.
"""
function square_free_mod_p(f::P, prime::Integer) where {C, D, P <: Polynomial{C, D}}
    f_poly = convert_polynomial(Val(prime), f)
    sqr_free = square_free(ZModP{C, prime}, f_poly)
    sqr_free = convert_polynomial(Val(prime), sqr_free)
    return sqr_free
end

#############################
# Queries about polynomials #
#############################

"""
Check if two polynomials are the same.
"""
function ==(p1::P, p2::P)::Bool where {C, D, P <: Polynomial{C, D}}
    if length(p1) != length(p2)
        return false
    end

    return all(t1 == t2 for (t1, t2) in zip(p1, p2))
end

"""
Check if a polynomial is equal to a constant `n`.
"""
function ==(p::Polynomial{C, D}, n::T) where {C, D, T <: Number}
    degree(p) != 0 && return false
    return leading(p).coeff == C(n)
end

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials (of the same concrete subtype).
"""
function -(p1::P, p2::P)::P where {C, D, P <: Polynomial{C, D}} 
    return p1 + (-p2)
end

"""
Multiplication of polynomial and term.
"""
function *(t::Term{C, D}, p::P)::P where {C, D, P <: Polynomial{C, D}} 
    return iszero(t) ? P() : P(map((pt)->t*pt, p))
end
*(p::P, t::Term{C, D}) where {C, D, P <: Polynomial{C, D}} = t*p

"""
Multiplication of polynomial and an integer.
"""
*(p::P, n::Integer) where {C, D, P <: Polynomial{C, D}} = Term{C, D}(C(n), D(0))*p
*(n::Integer, p::P) where {C, D, P <: Polynomial{C, D}} = p*n

"""
Integer division of a polynomial by an integer modulo a prime.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
function div_mod_p(p::P, n::Integer, prime::Integer) where {C, D, P <: Polynomial{C, D}}
    P( map((pt)->(div_mod_p(pt, n, prime)), p) )
end

"""
Take the mod of a polynomial with an integer.
"""
function mod(f::P, p::Int)::P where {C, D, P <: Polynomial{C, D}}
    f_out = P()
    for t in f
        new_t = mod(t, p)
        !iszero(new_t) && push!(f_out, new_t)
    end
    return trim!(f_out)
end

"""
Power of a polynomial mod prime.
"""
function pow_mod(p::P, n::Int, prime::Int) where {C, D, P <: Polynomial{C, D}}
    n < 0 && error("No negative power")
    
    # I do not know if this is cheating but this syntax and all the functions
    # make me dizzy so I implemented it in python and translated it

    # to return
    ans = one(p)

    # 0^n case
    if iszero(p)
        if n == 0
            return one(p)
        else
            return zero(p)
        end
    # 1^n
    elseif p == one(p)
        return one(p)
    
    # x^n case
    elseif p == x_poly(p)
        return P([Term{C, D}(one(C), D(n))])
    else
        base = p
        while n > 0
            # use a mask. if the bit is set we want to compute base^{(b_i)*2^j}
            # where b_i ∈ {0, 1} and j is the currently index bit of n
            if (n & 1) == 1
                # if the bit was set, multiply and mod it. If the bit was not set then
                # the power of the base has no contribution to the binary expansion
                # of the exp
                ans = mod(ans * base, prime)
            end

            # no matter what, divide by 2 since we want to look at each bit
            n >>= 1

            # repeated squaring with mod to control size
            base = mod(base * base, prime)
            
        end
    end
    return ans
end