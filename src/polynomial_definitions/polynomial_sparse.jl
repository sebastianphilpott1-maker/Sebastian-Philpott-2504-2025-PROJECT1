#############################################################################
#############################################################################
#
# This file defines the `PolynomialSparse` type with several operations 
#                                                                               
#############################################################################
#############################################################################

"""
The idea is simple, but I am a bit dumb so it took a while

Pass this function a vector of terms - however disorganised - and it
will remove zero terms, as well as walk over the vector, adding polynomials
of equivalent degree, and discarding them if they cancel to zero.

vec = [Term(0, 0), Term(2, 1), Term(0, 2), Term(2, 1)]
build_sparse_heap(vec) = Heap(Term(4, 1))

Aside: I am kind of convinced as of 4/10/2025 that this is wrong but things pass

We only store terms of unique degree. So if we see Term(2, 2) and Term(-2, 2)
we don't to cancel this out, right? Just store -2x. Or don't we?

Not sure but right now it works.
"""
function build_heap_sparse(vt::Vector{Term{C, D}}) where {C, D}
    # Hashmap so we can store/access potential terms of our sparse polynomial
    term_map = Dict{D, C}()

    # look at all the terms in the supplied list
    for term in vt

        # if its zero, move on
        iszero(term) && continue

        # store the coefficient and degree
        cur_coeff = term.coeff
        cur_degree = term.degree

        # we want to add constituent polynomial parts. If we try to get
        # something that doesn't exist yet, the get call defaults to zero
        # plus our new cur_coeff. we protected vs zero terms already
        cur_term_sum = get(term_map, cur_degree, zero(C)) + cur_coeff

        # if what we just calculated summed to zero (ie 3x - 3x) then delete
        # note that even if we have 3x - 3x + 3x in a row, we still delete
        # but is retrieved on the third loop
        if iszero(cur_term_sum)
            delete!(term_map, cur_degree)
        else
            term_map[cur_degree] = cur_term_sum
        end
    end

    # From the built map, assemble a new vector of terms.
    ret_vec_terms = Vector{Term{C, D}}()
    for (deg, coefficient) in term_map
        push!(ret_vec_terms, Term{C, D}(coefficient, deg))
    end

    # use Heap! to get a valid max heap
    return Heap!(ret_vec_terms)
end

#########################################
# PolynomialSparse type and construction #
#########################################

"""
A Polynomial type - designed to be for polynomials with integer coefficients.

"""
struct PolynomialSparse{C, D} <: Polynomial{C, D}
    terms::Heap{Term{C, D}}
    
    # zero sparse polynomial is an empty heap where the type passed to heap
    # is a term which will be stored in the heaps underlying vector/list
    PolynomialSparse{C, D}() where {C, D} = new(Heap(Term{C, D}))

    # Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparse{C, D}(vt::Vector{Term{C, D}}) where {C, D}

        # we use our helper function to construct a sparsepolynomial
        heapified_heap = build_heap_sparse(vt)
        new(heapified_heap)
    end

    # another inner constructor that lets us avoid building the dict map again
    # if we already know heap property is maintained
    function PolynomialSparse{C, D}(heap::Heap{Term{C, D}}) where {C, D}
        new(heap)
    end
end

##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.

Note that I don't think this works the way I think..? But added .data since the underlying
structure of the sparse polynomial is a heap whose underlying structure is a vector.
"""
function iterate(p::PolynomialSparse{C, D}, state = 1) where {C, D}
    return iterate(p.terms.data, state)
end

##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.

length automatically detects the length function in heap
"""
function length(p::PolynomialSparse{C, D}) where {C, D}
    return length(p.terms)
end

"""
The leading term of the polynomial.

Max heap implies biggest is at the top
"""
function leading(p::PolynomialSparse{C, D})::Term{C, D} where {C, D}
    isempty(p.terms) ? zero(Term{C, D}) : peek(p.terms) 
end

"""
The term of smallest degree in this polynomial.

Same logic, just do it over the array now. 

Unless I am mistaken, we cannot easily get the term of smallest degree
since we can't just donkey search from the back. How do we know when to stop?
"""
function last(p::PolynomialSparse{C, D}) where {C, D} 
    iszero(p) && return leading(p) # zero Term
    cur_min = p.terms.data[1]
    for term in p.terms.data
        if term.degree < cur_min.degree
            cur_min = term
        end
    end
    return cur_min
end


################################
# Pushing and popping of terms #
################################

"""
add func coz push got mucked up
"""
function add!(p::PolynomialSparse{C, D}, t::Term{C, D}) where {C, D}
    push!(p.terms, t)
    return p
end

"""
Push a new leading term into the polynomial (note - a constant can be pushed onto the zero polynomial).
"""
function push!(p::PolynomialSparse{C, D}, t::Term{C, D}) where {C, D}
    push!(p.terms, t)
    return p
end

"""
Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.

We rely on the heap implementation of pop
"""
function pop!(p::PolynomialSparse{C, D})::Term{C, D} where {C, D}
    iszero(p) && return zero(Term{C, D})
    return pop!(p.terms)
end

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

# We again can rely on the Julia type system - the implementations for `Polynomial` will work for `PolynomialSparse`

#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same.

Note - even though this is done for `Polynomial`, we can override it for `PolynomialSparse`
to leverage Julia's speed with vectors.

I can't remember what I had before but my equality absolutely didn't work so i used maps
"""
function ==(p1::PolynomialSparse{C, D}, p2::PolynomialSparse{C, D}) where {C, D}
    if length(p1) != length(p2)
        return false
    end
    # IM OVER IT. USE MAPS
    # https://www.youtube.com/watch?v=kVgy1GSDHG8
    trusty_hash_map = Dict{D, C}()

    # for every term in p1, map its degree to a coeff
    for t1 in p1
        trusty_hash_map[t1.degree] = t1.coeff
    end

    for t2 in p2
        trusty_variable = get(trusty_hash_map, t2.degree, nothing)

        # if there ever exists an equality check where the degree of t2
        # does not exist in the map built from terms of p1, we do not 
        # have equality. 
        if trusty_variable === nothing || trusty_variable != t2.coeff
            return false
        end
    end
    
    # surely this is efficient and works clueless
    return true
end

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

# Again, technically the abstract implementations for `Polynomial` will work correctly for `PolynomialSparse`.
# However, this doesn't mean they are particularly efficient - if you wish you can override any particular
# operation and re-implement it for `PolynomialSparse`.`