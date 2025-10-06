"""
Test product of polynomials.
"""
function pow_test_poly(::Type{P};
    N::Int = 10^2, power::Int = 11, seed::Int = 0
    ) where {C, D, P <: Polynomial{C, D}}
    Random.seed!(seed)
    
    for blah in 1:N
        p = rand(P; degree = 6, max_coeff = 6)
        ret = p^power
        
        # begin dummy power
        ans = leading(p)
        # start from 2 coz ans already first term 
        for _ in 2:power
            ans *= leading(p)
        end
        @assert leading(ret) == ans
        if (mod(blah, 50)) == 0
            println("Leading term of p^11 is ", leading(ret), " which is equal to leading term of p multiplied 11 times ", ans)
        end
    end

    for blah in 1:N
        p = rand(P; degree = 6, max_coeff = 6)
        ret = p^power

        # more dummy power
        ans = one(p)
        for _ in 1:power
            ans *= p
        end

        @assert ret == ans
        if (mod(blah, 50)) == 0
            println("If we multiply p^11 out it is equal to multiplying p by itself 11 times")
            println(ret == ans)
        end
    end

    println("pow_test_poly for $(P) - PASSED")
end

function pow_mod_test_poly(::Type{P};
    N::Int = 10^2, power::Int = 11, prime:: Int = 17, seed::Int = 0
    ) where {C, D, P <: Polynomial{C, D}}
    Random.seed!(seed)

    for blah in 1:N
        p = rand(P; degree = 6, max_coeff = 6)
        ret = pow_mod(p, power, prime)
        
        # begin dummy power
        ans = one(p)
        base = mod(p, prime)
        dummy_power = power
        while dummy_power > 0
            if (dummy_power & 1) == 1
                ans = mod(ans * base, prime)
            end
            dummy_power >>= 1
            base = mod(base * base, prime)
        end
        @assert ret == ans
        if (mod(blah, 50)) == 0
            println(leading(ret), " is the leading term of p to the power of 11 mod 17 is which is equal to ", leading(ans))
        end
    end

    for _ in 1:N
        p = rand(P; degree = 6, max_coeff = 6)
        ret = pow_mod(p, power, prime)

        # more dummy power
        ans = one(p)
        for blah in 1:power
            ans *= mod(p, prime)

            # mod the ans to control size
            ans = mod(ans, prime)
        end

        @assert ret == ans
    end

    println("passped for $(P)")
end