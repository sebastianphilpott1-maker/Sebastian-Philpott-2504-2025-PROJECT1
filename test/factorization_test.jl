#############################################################################
#############################################################################
#
# This file contains units tests for polynomial factorization
#                                                                               
#############################################################################
#############################################################################


"""
Test factorization of polynomials.
"""
function factor_mod_p_test_poly(::Type{P};
    N::Int = 10, seed::Int = 0, primes::Vector{Int} = [5, 7, 11]
    ) where {C, D, P <: Polynomial{C, D}}
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")

            # Ensure our random polynomial does not have vanishing leading term mod prime
            p = rand(P; monic = true)

            factorization = factor_mod_p(p, prime)
            pr = expand_factorization(factorization)
            @assert mod(p-pr, prime) == 0
        end
    end

    println("\n\nfactor_mod_p_test_poly for $(P) - PASSED")
end

