#############################################################################
#############################################################################
#
# This file implements polynomial addition for sparse polynomials.
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.

If we want to add a term in the sparse polynomial, we need to check
for its existence.

Honestly at this point I don't even know if thats true.
Sometimes I feel like we can add a term with lower degree than leading
and other times I think we can't.
"""
function +(p::PolynomialSparse{C, D}, t::Term{C, D}) where {C, D}
    iszero(t) && return p
    iszero(p) && return PolynomialSparse{C, D}([t])
    data = p.terms.data

    # A flag to tell us if we didn't find the degree in the polynomial
    term_doesnt_exist = true

    # start looking for term of matching degree
    for i in 1:length(p)
        if data[i].degree == t.degree

            # these are terms, add them
            potential_term = data[i] + t

            # Check if it is zero
            if iszero(potential_term)

                # If the potential updated term is zero 
                # we need to remove it from the polynomial
                # and restore heap ordering
                delete_at!(p.terms, i)
            
            # If it isn't zero, overide the existing element
            # No maintenance of heap order necessary since
            # we order via degree
            else
                data[i] = potential_term
            end
            term_doesnt_exist = false
            # after this check we break since we don't want to keep
            # scanning since there are no dupes
            break
        end
    end

    # if a term of equal degree doesn't exist, add the desired term
    if term_doesnt_exist
        add!(p, t)
    end
    return p
end

function +(p1::PolynomialSparse{C, D}, p2::PolynomialSparse{C, D})::PolynomialSparse{C, D} where {C, D}
    # This is probably causing overflow issues
    p = deepcopy(p1)
    q = deepcopy(p2)

    ans = PolynomialSparse{C, D}()


    while !iszero(p) || !iszero(q)
        if degree(p) == degree(q)
            ltp, ltq = pop!(p), pop!(q)
            if !iszero(ltp + ltq)
                add!(ans, ltp + ltq)
            end
        
        elseif degree(p) > degree(q)
            add!(ans, pop!(p))

        else
            add!(ans, pop!(q))
        end
    end
    return ans
end


+(p::PolynomialSparse{C, D}, n::C) where {C, D} = p + Term{C, D}(n, 0)
+(n::C, p::PolynomialSparse{C, D}) where {C, D} = Term{C, D}(n, 0) + p