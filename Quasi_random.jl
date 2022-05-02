"""
    Harmonious_number(n::Int)

Harmonious number (or clarity ratio) generator
for n-dimensional Euclidean space.

For given natural number 2 ≥ n ,let 

    fₙ(x) = xⁿ - x - 1.

Then there exists real number χₙ∈(1,2) such that fₙ(χₙ)=0 
and for every root r ∈ C of the polynomial fₙ the following
statement holds:

    r ≠ χₙ ⇨ |r| < χₙ

For given natural number 2 ≥ n ,let

    Hₖ⁽ⁿ⁾ = 1,                        k = 1,2,...,n,
            Hₖ₋ₙ⁽ⁿ⁾ + Hₖ₋ₙ₊₁⁽ⁿ⁾,      k>n 

Then nuber χₙ is defined with

    χₙ = limₖ→∞ (Hₖ₊₁⁽ⁿ⁾ / Hₖ⁽ⁿ⁾).

Marohnić, L. and Strmečki, T., 2016, November. 
Plastic number: construction and applications. 
In 5th Virtual International Conference on Advanced Research 
in Scientific Areas (ARSA-2016) Slovakia.

Returns:
    χₙ: harmonious number given n

""" 
function Harmonious_number(n::Int)

    n = n+1
    N = n * 500 #number of sequence
    Hₙ = ones(N+1)
    k = n+1

    while k <= N+1

        Hₙ[k] = Hₙ[k-n] + Hₙ[k-n+1]
        k +=1
    end

    χₙ = Hₙ[N+1] / Hₙ[N]

    return χₙ
end

"""
    function R_n_sequence(s₀::Real, d::Int, N::Int)

recurrence R-sequence that falls in the category of those 
sequences that are based on specially chosen irrational numbers. 

The following parameter-free d-dimensional open(infinite) sequence, R_d(ϕ_d) 
has excellent low discrepancy characteristics when compared to other existing methods.

    tₙ = {s₀ + nα},         n=1,2,3,...
    α = (1/ϕ_d, 1/ϕ²_d, 1/ϕ³_d, ..., 1/ϕᵈ_d)
    ϕ is the Harmonic number 

http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/

"""
function R_n_sequence(s₀::Real, d::Int, N::Int)

    ϕ = Harmonious_number(d)
    α = zeros(d)
    t_n = zeros(Float64, (d, N))

    for i in 1:d
        α[i] = (1 / ϕ^i)
    end

    for j=1:N
        t_n[:, j] = (s₀ .+ j .* α) .% 1
    end
    
    return t_n
end