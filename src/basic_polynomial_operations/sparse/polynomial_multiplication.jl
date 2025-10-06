#############################################################################
#############################################################################
#
# This file implements polynomial multiplication for dense polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Multiplcation should preserve order i think
"""
function *(p::PolynomialSparse{C, D}, t::Term{C, D}) where {C, D}
    # term x term multiplication is defined
    # scalar multiplication preserves order
    # From this, use convenience constructor that assumes ordering so 
    # we don't have to build a map again
    return PolynomialSparse{C, D}(map_heap(p.terms, term -> term * t))
end

"""
Going back over and changing all the ops for n::Int to n::C

If i understand, since the ZModP object is a integer subtype, then
scalar multiplciation like this should be defined? I think?
"""
*(p::PolynomialSparse{C, D}, n::Integer) where {C, D} = p * Term(n, 0)
*(n::Integer, p::PolynomialSparse{C, D}) where {C, D} = Term(n, 0) * p

"""
Will comment code but idea is this:

In general, product of two polynomials p1 and p2 with terms xi, yj is

x1(y1 + y2 + ...) + x2(y1 + y2 + ...) etc

All of these are terms, for which we have defined multiplciation for. Moreover,
we already have addition defined so the result is as follows.

Ended up using a hashmap implementation. Because... I actually don't know anymore
I think I had errors or I had overflow issues with using + between polynomials
because I called deepcopy too many times.
"""
function *(p1::PolynomialSparse{C, D}, p2::PolynomialSparse{C, D}) where {C, D}
    # "Double for loop" over the terms to do it cross-wise
    # All we are doing is getting a term in p1, adding to answer the product
    # of x1 and y1, then add prod of x1 and y2 etc etc
    # but now that I think about it its slow as dog because we loop in + so
    # keeping this here because it "works" but it deepcopies 9 billion things
    # for t in p1_terms
    #     new_summand = p2 * t
    #     # our addition function will handle everything ::Clueless::
    #     ans = ans + new_summand
    # end
    # return ans

    # throw a hashmap at it - use deg, coeff pair since degree is unique
    mult_mapping = Dict{D, C}()

    # for every term in p1, look at every term in p2
    for t1 in p1
        for t2 in p2

            # store a new terms potential degree and coeff
            cur_deg = t1.degree + t2.degree
            cur_coeff = t1.coeff * t2.coeff

            # if the sum is zero, don't add it and move on
            if iszero(cur_coeff)
                continue
            end

            # If coeff != 0, at the degree key, sum whatever exists with the update
            # If location doesn't exist, default get 0 and add cur_coeff
            mult_mapping[cur_deg] = get(mult_mapping, cur_deg, zero(C)) + cur_coeff
        end
    end

    new_terms = Term{C, D}[]
    # loop over key value pairs
    for (degree, coeff) in mult_mapping

        # if a zero coeff slipped in, do not add it.
        if iszero(coeff)
            continue
        else
            push!(new_terms, Term{C, D}(coeff, degree))
        end
    end
    return PolynomialSparse{C, D}(new_terms)
end

# function ^(p::PolynomialSparse{C, D}, n::Integer) where {C, D}
#     """
#     Ripping my code from elsewhere
#     """
#     # If the power is negative, error
#     n < 0 && error("No negative power")

#     # Empty constructor to store our answer
#     ans = one(p)

#     # 0^n case
#     if p == zero(p)
#         if n == 0
#             return one(p)
#         else
#             return zero(p)
#         end
#     # 1^n
#     elseif p == one(p)
#         return one(p)
    
#     # x^n case
#     elseif p == x_poly(p)
#         return PolynomialSparse{C, D}([Term{C, D}(one(C), D(n))])
#     else
#         base = p
#         while n > 0
#             if (n & 1) == 1
#                 ans = ans * base
#             end
#             n >>= 1
#             base = base * base
#         end
#     end
#     return ans
# end



