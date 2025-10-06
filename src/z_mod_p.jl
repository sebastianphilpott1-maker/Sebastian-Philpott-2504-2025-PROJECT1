using Primes
"""
For my own sanity, we are implementing Z_{p}, which is the set of integers
mod p.

Ex 1) Z_{5} = {0, 1, 2, 3, 4}
Ex 2) Z_{9} = {0, 1, 2, 3, 4, 5, 6, 7, 8}

Our goal - to my understanding - is to model the elements of Z_{p} where
p is a prime number.

Formally, each time a ZModP object is constructed, yes, it is just an integer
however it represents the equivalance class of integers congruent to val mod N

"""
struct ZModP{T <: Integer, N} <: Integer
    val::T

    # val ∈ [0,...,N - 1]
    # T is the type of val, N is the prime p in Zp

    # Inner constructor that enforces the invariant that 0 ≤ val ≤ N - 1
    # which gives us a number in [0, N-1]
    function ZModP{T, N}(a::T) where {T <: Integer, N}
        !isprime(N) && error("N must be prime!")
        return new(mod(a, N))
    end
end

"""
Ran into issues when testing of this variety:

ERROR: MethodError: no method matching ZModP{BigInt, 5}(::Int64)
The type `ZModP{BigInt, 5}` exists, but no method is defined for 
this combination of argument types when trying to construct it.

Now, this convenience constructor reads any input and just donkey
converts it to type T
"""
function ZModP{T, N}(a::Integer) where {T <: Integer, N}
    return ZModP{T, N}(T(a))
end

"""
Would have been handy if i committed earlier but i was using this
every where and it was kind of pointless so i redid all OPS
without it but it lurks in the code
"""
function get_zmod_val(a::ZModP{T, N})::T where {T <: Integer, N}
    return a.val
end

function get_zmod_mod(::Type{ZModP{T, N}}) where {T, N}
    return N
end

modulus_of(::Type{ZModP{T,N}}) where {T,N} = N

function ZModP(a::T, ::Val{N})::ZModP{T, N} where {T <: Integer, N}
    return ZModP{T, N}(a)
end

function zero(::Type{ZModP{T, N}})::ZModP{T, N} where {T <: Integer, N}
    return ZModP{T, N}(zero(T))
end

function one(::Type{ZModP{T, N}})::ZModP{T, N} where {T <: Integer, N}
    return ZModP{T, N}(one(T))
end

# This took me a stupid long time to figure out bless the docs
# https://docs.julialang.org/en/v1/base/numbers/#Base.zero
function zero(a::ZModP{T, N})::ZModP{T, N} where {T <: Integer, N}
    return zero(typeof(a))
end

function one(a::ZModP{T, N})::ZModP{T, N} where {T <: Integer, N}
    return one(typeof(a))
end

function show(io::IO, a::ZModP{T, N}) where {T <: Integer, N}
    print(io, a.val)
end

function ==(a::ZModP{T, N}, b::ZModP{T, N})::Bool where {T <: Integer, N}
    return a.val == b.val
end

function ==(a::ZModP{T, N}, b::S)::Bool where {T <: Integer, S <: Integer, N}
    return a.val == b
end

function ==(a::S, b::ZModP{T, N})::Bool where {T <: Integer, S <: Integer, N}
    return b.val == a
end

####################
# BASIC STUFF #
####################

"""
Since Z_{p} is a field, we know it is closed under addition. Hence, we have
x = ZModP{BigInt, 5}(BigInt(11))
y = ZModP{BigInt, 5}(BigInt(13))

⟹ x + y = (11 + 13) % 5

ie. (11 mod 5 + 13 mod 5) mod 5 = (11 + 13) mod 5
"""
function +(a::ZModP{T, N}, b::ZModP{T, N}) where {T <: Integer, N}
    return ZModP{T, N}(a.val + b.val)
end

"""
Addition in a field Z_p between a ZModP object and integer
"""
function +(a::ZModP{T, N}, b::S) where {T <: Integer, S <: Integer, N}
    # use addition defined for two zmodp objs
    return a + ZModP{T, N}(T(b))
end

"""
Addition in a field Z_p between an integer and ZModP object

Not sure why we would implement it as described in the spec so until
I see an error, I am doing it like above
"""
function +(a::S, b::ZModP{T, N}) where {T <: Integer, S <: Integer, N}
    return ZModP{T, N}(T(a)) + b
end

"""
Take the negative (additive inverse in Z_{p}) of a ZModP object
"""
function -(a::ZModP{T, N}) where {T <: Integer, N}
    return ZModP{T, N}(-a.val)
end

"""
Modular subtraction between two Z_p objects
"""
function -(a::ZModP{T, N}, b::ZModP{T, N}) where {T <: Integer, N}
    return ZModP{T, N}(a.val - b.val)
end

"""

"""
function -(a::ZModP{T, N}, b::S) where {T <: Integer, S <: Integer, N}
    return a - ZModP{T, N}(T(b))
end


"""

"""
function -(a::S, b::ZModP{T, N}) where {T <: Integer, S <: Integer, N}
    return ZModP{T, N}(T(a)) - b
end

"""
"""
function *(a::ZModP{T, N}, b::ZModP{T, N}) where {T <: Integer, N}
    return ZModP{T, N}(a.val * b.val)
end

function *(a::ZModP{T, N}, b::S) where {T <: Integer, S <: Integer, N}
    return a * ZModP{T, N}(T(b))
end

function *(a::S, b::ZModP{T, N}) where {T <: Integer, S <: Integer, N}
    return ZModP{T, N}(T(a)) * b
end

"""
This seems right
"""
function inv(a::ZModP{T, N})::ZModP{T, N} where {T <: Integer, N}
    iszero(a) && throw(DivideError())
    return ZModP{T, N}(int_inverse_mod(T(a.val), T(N)))
end

function ÷(a::ZModP{T, N}, b::ZModP{T, N})::ZModP{T, N} where {T <: Integer, N}
    iszero(b) && throw(DivideError())
    return a * inv(b)
end

function ÷(a::ZModP{T, N}, b::S)::ZModP{T, N} where {T <: Integer, S <: Integer, N}
    # we need to ensure we return a ZModP and enforce the invariant upon the
    # result. We rely on our other division method since it uses inv already
    return a ÷ ZModP{T, N}(T(b))
end

function ÷(a::S, b::ZModP{T, N})::ZModP{T, N} where {T <: Integer, S <: Integer, N}
    return ZModP{T, N}(T(a)) ÷ b
end

"""
This function was written like this because I thought we had to compute
negative powers.

However, https://edstem.org/au/courses/24848/discussion/2956967 seemed to 
suggest otherwise but in any case inv has been imported to the project
"""
function ^(a::ZModP{T, N}, n::S) where {T <: Integer, S <: Integer, N}
    # similar to powmod handle cases
    if iszero(a)
        n > 0 && return zero(a)
        n == 0 && return one(a)
        n < 0 && error("negative exponent")
    end
    
    # guard against negative exponents
    n == 0 && return one(a)

    # Ensure it is in the range
    idx_copy = mod(n, N - 1)
    base = a
    ans = one(a)
    while idx_copy > 0
        if (idx_copy & 1) == 1
            ans = ans * base
        end
        idx_copy >>= 1
        base = base * base
    end
    return ans
end

"""
IMPLEMENT LATER THIS IS RUINING MY HACKY PRINTING IMPLEMENTATION
"""
# function abs(a::ZModP{T, N})::ZModP{T, N} where {T <: Integer, N}
#     return get_zmod_val(a)
# end

"""
I don't know why I've got any of these but they dont't really work
or I don't know why use them
"""
BigInt(a::ZModP{T, N}) where {T <: Integer, N} = BigInt(get_zmod_val(a))
Int(a::ZModP{T, N}) where {T <: Integer, N} = Int(get_zmod_val(a))

function get_zmodp_types(a::ZModP{T, N}) where {T <: Integer, N}
    return T, N
end