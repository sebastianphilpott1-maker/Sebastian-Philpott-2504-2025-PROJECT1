#############################################################################
#############################################################################
#
# This file implements polynomial multiplication for abstract polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Multiply two polynomials (of the same concrete subtype).
"""
function *(p1::P, p2::P)::P where {C, D, P <: Polynomial{C, D}}
    p_out = P()
    for t in p1
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end

"""
Power of a polynomial.

Updated method using repeated squaring, not controlling size with mod obvs
"""
function ^(p::P, n::Integer)::P where {C, D, P <: Polynomial{C, D}}
    n < 0 && error("No negative power")
    ans = one(p)
    # 0^n case
    if p == zero(p)
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
            if (n & 1) == 1
                ans = ans * base
            end
            n >>= 1
            base = base * base
        end
    end
    return ans
end

