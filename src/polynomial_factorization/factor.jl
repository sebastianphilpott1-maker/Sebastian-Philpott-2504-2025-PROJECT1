#############################################################################
#############################################################################
#
# This file implements factorization 
#                                                                               
#############################################################################
#############################################################################

"""
Factors a polynomial over the field Z_p.

Returns a vector of tuples of (irreducible polynomials (mod p), multiplicity) such that their product of the list (mod p) is f. 
Irreducibles are fixed points on the function `factor`.
"""
function factor_mod_p(f::P, prime::Integer)::Vector{Tuple{P, Integer}} where {C, D, P <: Polynomial{C, D}}
    if C <: ZModP
        return factor(C, f)
    else
        f_poly = convert_polynomial(Val(prime), f)
        zmodp_ret = factor(ZModP{C, prime}, f_poly)
        
        converted_ret = Tuple{P, Integer}[]
        
        # convert back - val param doesn't matter
        for i in 1:length(zmodp_ret)
            push!(converted_ret, (convert_polynomial(Val(prime), zmodp_ret[i][1]), zmodp_ret[i][2]))
        end
        return converted_ret
    end
end

"""
Expand a factorization.
"""
function expand_factorization(factorization::Vector{Tuple{P, Integer}})::P where {P <: Polynomial}
    length(factorization) == 1 && return first(factorization[1])^last(factorization[1])
    return *([first(tt)^last(tt) for tt in factorization]...)
end

"""
Compute the number of times g divides f modulo a prime.
"""
function multiplicity_mod_p(f::P, g::P, prime::Integer)::Integer where {C, D, P <: Polynomial{C, D}}
    if C <: ZModP
        return multiplicity(C, f, g)
    else
        f_poly = convert_polynomial(Val(prime), f)
        g_poly = convert_polynomial(Val(prime), g)
        return multiplicity(ZModP{C, prime}, f_poly, g_poly)
    end
end


"""
Distinct degree factorization.

Given a square free polynomial `f` returns a list, `g` such that `g[k]` is a product of irreducible 
polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2, such that the product of the list (mod `prime`) 
is equal to `f` (mod `prime`).
"""
function dd_factor_mod_p(f::P, prime::Integer)::Array{P} where {C, D, P <: Polynomial{C, D}}
    if C <: ZModP
        return dd_factor(C, f)
    else
        f_poly = convert_polynomial(Val(prime), f)
        zmodp_ret = dd_factor(ZModP{C, prime}, f_poly)
        converted_ret = Array{P}(undef, length(zmodp_ret))

        # convert back - val param doesn't matter
        for i in 1:length(zmodp_ret)
            push!(converted_ret, convert_polynomial(Val(prime), zmodp_ret[i]))
        end
        return converted_ret
    end
end

"""
Distinct degree split.

Returns a list of irreducible polynomials of degree `d` so that the product of that list (mod prime) is the polynomial `f`.
"""
function dd_split_mod_p(f::P, d::Integer, prime::Integer)::Vector{P} where {C, D, P <: Polynomial{C, D}}
    if C <: ZModP
        return dd_split(C, f, d)
    else
        f_poly = convert_polynomial(Val(prime), f)
        zmodp_ret = dd_split(ZModP{C, prime}, f_poly, d)
        converted_ret = Vector{P}()
        for i in 1:length(zmodp_ret)
            push!(converted_ret, convert_polynomial(Val(prime), zmodp_ret[i]))
        end
        return converted_ret
    end
end