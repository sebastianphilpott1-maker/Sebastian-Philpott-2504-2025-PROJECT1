#############################################################################
#############################################################################
#
# This file implements polynomial addition for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.

This must be overridden by concrete subtypes.
"""
function +(p::Polynomial{C, D}, t::Term{C, D}) where {C, D}
    not_implemented_error(p, "Polynomial + Term")
end
+(t::Term{C, D}, p::Polynomial{C, D}) where {C, D} = p + t

"""
Add two polynomials of the same concrete subtype.

Note: This operation may be slow for some concrete subtypes. You may wish to override this to factor
in the details of your polynomial representation when implementing your concrete subtype.
"""
function +(p1::P, p2::P)::P where {C, D, P <: Polynomial{C, D}}
    p = deepcopy(p1)
    for t in p2
        p += t
    end
    return p
end

"""
Add a polynomial and an integer.
"""
+(p::Polynomial{C, D}, n::Integer) where {C, D} = p + Term{C, D}(C(n), D(0))
+(n::Integer, p::Polynomial{C, D}) where {C, D} = p + Term{C, D}(C(n), D(0))