#############################################################################
#############################################################################
#
# This file defines the Term type with several operations 
#                                                                               
#############################################################################
#############################################################################

##################################
# Term type and utility function #
##################################

"""
A term.

For Term{C, D},
    C: The type of the coefficient of the term.
    D: The type of the degree of the term.
"""
struct Term{C <: Integer, D <: Integer} #structs are immutable by default
    coeff::C
    degree::D

    function Term{C, D}(coeff::C, degree::D) where {C, D}
        degree < 0 && error("Degree must be non-negative")
        new(coeff, degree)
    end
end

######################
# Outer Constructors #
######################

""" Convenience outer constructor to infer coefficient/degree types. """
function Term(coeff::C, degree::D) where {C, D}
    Term{C, D}(coeff, degree)
end


"""
Creates the zero term.
"""
function zero(::Type{Term{C, D}})::Term{C, D} where {C, D} 
    Term(zero(C), zero(D))
end
zero(t::Term) = zero(typeof(t))

"""
Creates the unit term.
"""
function one(::Type{Term{C, D}})::Term{C, D} where {C, D}
    Term(one(C), zero(D))
end
one(t::Term) = one(typeof(t))

###########
# Display #
###########

"""
Show a term.
"""
show(io::IO, t::Term) = print(io, "$(t.coeff)x^$(t.degree)")

########################
# Queries about a term #
########################

"""
Check if a term is 0.
"""
iszero(t::Term)::Bool = iszero(t.coeff)

"""
Compare two terms (with the same coefficient/degree types).
"""
function isless(t1::Term{C, D}, t2::Term{C, D})::Bool  where {C, D}
    t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  
end


function (==)(t1::Term{C, D}, t2::Term{C, D}) where {C, D}
    # Without this zmodp comparison fails
    if iszero(t1) && iszero(t2)
        return true
    end
    t1.coeff == t2.coeff && t1.degree == t2.degree
end

"""
Evaluate a term at a point x.
"""
function evaluate(t::Term, x::S) where {S <: Number} 
    t.coeff * x^t.degree
end

##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree (with the same coefficient/degree types).
"""
function +(t1::Term{C, D}, t2::Term{C, D})::Term{C, D} where {C, D}
    @assert t1.degree == t2.degree
    Term{C, D}(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
function -(t::Term{C, D},)::Term{C, D} where {C, D} 
    Term{C, D}(-t.coeff, t.degree)  
end

"""
Subtract two terms with the same degree (and the same coefficient/degree types).
"""
function -(t1::Term{C, D}, t2::Term{C, D})::Term{C, D} where {C, D}
    t1 + (-t2) 
end

"""
Multiply two terms (with the same coefficient/degree types).
"""
function *(t1::Term{C, D}, t2::Term{C, D})::Term{C, D} where {C, D}
    Term{C, D}(t1.coeff * t2.coeff, t1.degree + t2.degree)
end

"""
Multiply a term by a constant.
"""
function *(t::Term{C, D}, n::S)::Term{C, D} where {C, D, S <: Integer}
    Term{C, D}(t.coeff * n, t.degree)
end
function *(n::S, t::Term{C, D})::Term{C, D} where {C, D, S <: Integer}
    t * n
end

"""
Power of a term.
"""
function ^(t::Term{C, D}, n::S)::Term{C, D} where {C, D, S <: Integer}
    Term{C, D}(t.coeff^n, t.degree*D(n))
end

"""
Compute the symmetric mod of a term with an integer.
"""
function mod(t::Term{C, D}, p::Integer)::Term{C, D} where {C <: Integer, D}
    if C <: ZModP
        return t 
    else
        Term{C, D}(mod(t.coeff, p), t.degree)
    end
end

"""
Compute the derivative of a term.
"""
function derivative(t::Term{C, D})::Term{C, D} where {C, D}
    Term{C, D}(t.coeff * C(t.degree), max(t.degree - one(D), zero(D)))
end

"""
Exact division when the coefficient type `C` is a field (e.g. Zp).

You will need to override this in Task 5 for ZModP specifically.
Do NOT modify this particular version of the function.

Since division is defined for our zmodp object and it is indeed a field
we can just call our defined function and julia will do the rest for us

Initially wrong: Was just doing division and return numbers, not terms
"""
function div(t1::Term{C, D}, t2::Term{C, D}) where {C, D}
    @assert t1.degree >= t2.degree
    if iszero(t1.coeff ÷ t2.coeff)
        return zero(Term{C, D})
    end
    return Term{C, D}(t1.coeff ÷ t2.coeff, D(t1.degree - t2.degree))
end

"""
Divide two terms (with the same coefficient/degree types). 
Returns a function of an integer.

Note: You will need to override this for division where the coefficients are of type ZModP (Task 5)
There we can do exact division, so we can simply do `t1.coeff ÷ t2.coeff`.
"""
function div_mod_p(t1::Term{C, D}, t2::Term{C, D}, prime::Integer) where {C, D}
    @assert t1.degree ≥ t2.degree
    # convert from input type Term{C, D} to Term{ZModP{C, N}, D}
    # anti dizzy comment: term constructor -> Term{C, D} where 
    # C -> ZModP{C, prime}
    # Term param (C, D) -> (ZModP{C, prime}, D)
    # ZModP{C, Prime} also needs param which is coeff of term

    # ---TERM CTOR -------TERM PARAM/ZMOD CTOR---ZMOD PARAM--D PART OF TERM
    t1_z = Term{ZModP{C, prime}, D}(ZModP{C, prime}(t1.coeff), t1.degree)
    t2_z = Term{ZModP{C, prime}, D}(ZModP{C, prime}(t2.coeff), t2.degree)
    
    # call div
    zmod_obj = div(t1_z, t2_z)

    # coeff is the value of the resultant zmod object
    return Term{C, D}(C(zmod_obj.coeff), D(t1.degree - t2.degree))
end

"""
Integer divide a term by an integer.
"""
function div_mod_p(t::Term{C, D}, n::Integer, prime::Integer) where {C <: Integer, D} 
    return div_mod_p(t, Term(C(n), zero(D)), prime)
end

#############################
# Vectorization with a term #
#############################

""" Enable broadcasting an operation on a term """
broadcastable(t::Term) = Ref(t)
