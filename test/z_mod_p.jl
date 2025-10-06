function test_zmodp_obj(N::Int = 1000)
    println("Begin comprehensive testing of zmodp")
    test_zmodp_add(N)
    test_zmodp_subtract(N)
    test_zmodp_mult(N)
    test_zmodp_division(N)
    test_zmodp_inv(N)
    test_zmodp_exp(N)
    println("PASSED ALL TESTS YAY")
end

function test_zmodp_add(N::Int)
    Random.seed!(0)
    primes = [2,3,5,7,11,13,17,19,23,29]

    for p in primes
        for blah in 1:N
            # create two ints in the range 10^-4, 10^4
            int_a = rand(-10^4:10^4)
            int_b = rand(-10^4:10^4)

            # make some zmodp objects to test
            zed_a = ZModP{BigInt, p}(BigInt(int_a))
            zed_b = ZModP{BigInt, p}(BigInt(int_b))

            @assert (zed_a + zed_b) == mod(int_a + int_b, p)
            @assert (zed_b + zed_a) == mod(int_b + int_a, p)
            @assert (int_a + zed_b) == mod(int_a + int_b, p)
            @assert (zed_a + int_b) == mod(int_a + int_b, p)

            if mod(blah, 500) == 0 && p == 23
                println("For prime = $(p), ", " $(zed_a) + $(zed_b) mod $(p) = ", (zed_a + zed_b), " which is the same as mod($(int_a) + $(int_b), p) = ", mod(int_a + int_b, p))
            end
        end
    end
    println("test_zmodp_add - PASSED")
end

function test_zmodp_subtract(N::Int)
    Random.seed!(0)
    primes = [2,3,5,7,11,13,17,19,23,29]

    for p in primes
        for blah in 1:N
            # create two ints in the range 10^-4, 10^4
            int_a = rand(-10^4:10^4)
            int_b = rand(-10^4:10^4)

            # make some zmodp objects to test
            zed_a = ZModP{BigInt, p}(BigInt(int_a))
            zed_b = ZModP{BigInt, p}(BigInt(int_b))

            @assert (zed_a - zed_b) == mod(int_a - int_b, p)
            @assert (zed_b - zed_a) == mod(int_b - int_a, p)
            @assert (int_a - zed_b) == mod(int_a - int_b, p)
            @assert (zed_a - int_b) == mod(int_a - int_b, p)
            @assert (-zed_a) == mod(-int_a, p)

            if mod(blah, 500) == 0 && p == 17
                println("For prime = $(p), ", " $(zed_a) - $(zed_b) mod $(p) = ", (zed_a - zed_b), " which is the same as mod($(int_a) - $(int_b), p) = ", mod(int_a - int_b, p))
            end
        end
    end
    println("test_zmodp_subtract - PASSED")
end

function test_zmodp_mult(N::Int)
    Random.seed!(0)
    primes = [2,3,5,7,11,13,17,19,23,29]

    for p in primes
        for blah in 1:N
            # create two ints in the range 10^-4, 10^4
            int_a = rand(-10^4:10^4)
            int_b = rand(-10^4:10^4)

            # make some zmodp objects to test
            zed_a = ZModP{BigInt, p}(BigInt(int_a))
            zed_b = ZModP{BigInt, p}(BigInt(int_b))

            @assert (zed_a * zed_b) == mod(int_a * int_b, p)
            @assert (zed_b * zed_a) == mod(int_b * int_a, p)
            @assert (int_a * zed_b) == mod(int_a * int_b, p)
            @assert (zed_a * int_b) == mod(int_a * int_b, p)

            if mod(blah, 500) == 0 && p == 13
                println("For prime = $(p), ", " $(zed_a) * $(zed_b) mod $(p) = ", (zed_a * zed_b), " which is the same as mod($(int_a) * $(int_b), p) = ", mod(int_a * int_b, p))
            end
        end
    end
    println("test_zmodp_mult - PASSED")
end

function test_zmodp_division(N::Int)
    Random.seed!(0)
    primes = [2,3,5,7,11,13,17,19,23,29]

    for p in primes
        for blah in 1:N
            # create two ints in the range 10^-4, 10^4
            int_a = rand(-10^4:10^4)
            int_b = rand(-10^4:10^4)

            # make some zmodp objects to test
            zed_a = ZModP{BigInt, p}(BigInt(int_a))
            zed_b = ZModP{BigInt, p}(BigInt(int_b))
            
            # guard against zero division
            if mod(int_b, p) != 0
                @assert (zed_a ÷ zed_b) == mod(int_a * int_inverse_mod(int_b, p), p)
                @assert (int_a ÷ zed_b) == mod(int_a * int_inverse_mod(int_b, p), p)
                @assert (zed_a ÷ int_b) == mod(int_a * int_inverse_mod(int_b, p), p)
            end
            if mod(int_a, p) != 0
                @assert (zed_b ÷ zed_a) == mod(int_b * int_inverse_mod(int_a, p), p)
            end
        end
    end
    println("test_zmodp_division - PASSED")
end

function test_zmodp_inv(N::Int)
    Random.seed!(0)
    primes = [2,3,5,7,11,13,17,19,23,29]

    for p in primes
        for blah in 1:N
            # create two ints in the range 10^-4, 10^4
            int_a = rand(-10^4:10^4)

            # make some zmodp objects to test
            zed_a = ZModP{BigInt, p}(BigInt(int_a))
            if mod(int_a, p) != 0
                @assert inv(zed_a) == mod(int_inverse_mod(int_a, p), p)
            end
        end
    end
    println("test_zmodp_inv - PASSED")
end

function test_zmodp_exp(N::Int)
    Random.seed!(0)
    primes = [2,3,5,7,11,13,17,19,23,29]

    for p in primes
        for blah in 1:N
            # create two ints in the range 10^-4, 10^4
            int_a = rand(-10^4:10^4)

            # make some zmodp objects to test
            zed_a = ZModP{BigInt, p}(BigInt(int_a))

            exponent = rand(0:10)
            # im not doing fast exp again cba
            # https://www.jlhub.com/julia/manual/en/function/powermod
            @assert (zed_a^exponent) == powermod(mod(int_a, p), exponent, p)
            if mod(blah, 500) == 0 && p == 29
                println("For prime = $(p), $(zed_a)^$(exponent) mod $(p) is ", zed_a^exponent, " which equals ", powermod(mod(int_a, p), exponent, p))
            end
        end
    end
    println("test_zmodp_exp - PASSED")
end

