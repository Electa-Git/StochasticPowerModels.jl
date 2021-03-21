"""
    APPF()

Accelerated Probabilistic Power Flow (APPF) is used to accelerate any sampling-
based PPF solver. As the APPF runs, it concurrently generates a low-dimensional
subspace of orthonormalized solution vectors. This subspace is used to construct
and update a reduced order model of the full nonlinear system, resulting in a 
highly efficient simulation for future voltage profiles.


```bibtex
Accelerated Probabilistic Power Flow via Model Order Reduction and Neumann 
Series Expansion, by S. Chevalier, L. Schenato, and L. Daniel, 2020.
```

Input: 
- S:  the sampled load points.
- x₀: initial voltage guess
- s₀: nominal power injection
- Jᶜ: reduced power flow Jacobian evaluated at x₀
- Hᶜ: reduced power flow Hessian evaluated at x₀
- 

"""
function APPF(L, U, P, S, x₀, s₀, Jᶜ, Hᶜ)
    xᶜ  = CartesianFromPolar().(x₀)
    x̂ᶜ  = norm(xᶜ)
    δx̂ᶜ = 0.0
    V   = xᶜ / x̂ᶜ
    Ĵ   = Jᶜ * V
    ŝ₀  = Ĵ' * s₀
    Ĝ   = Ĵ' * Ĵ
    Ĥ   = Ĵ' * Hᶜ * kron(V ,V)
    for Sⁱ in S
        Ŝ = Ĵ' * Sⁱ
        δx̂ᶜ, x = RMS(V, ŝ₀, Ĝ, Ĥ, x̂ᶜ, δx̂ᶜ, Ŝ)
        if norm(s(x) - S, Inf) > 1e-5
            NSBPF!(L, U, P, Sⁱ, x)
            DSE!(V, Ĵ, Ĝ, Ĥ, δx̂ᶜ, ŝ₀, s₀, x, Jᶜ, Hᶜ)
        end
    end 
    return x
end
function FEBS(L, U, P, b)
    # solve L * z = P' * b
    y = U \ z
    x = P * y
end
function RMS(V, ŝ₀, Ĝ, Ĥ, x̂ᶜ, δx̂ᶜ, Ŝ)
    ĝ = ŝ₀ + Ĝ * δx̂ᶜ + 1/2 * Ĥ * kron(δx̂ᶜ, δx̂ᶜ) - Ŝ
    while norm(ĝ, Inf) > 1e-5
        δx̂ᶜ .-= inv(Ĝ) * ĝ
        ĝ   =  ŝ₀ + Ĝ * δx̂ᶜ + 1/2 * Ĥ * kron(δx̂ᶜ, δx̂ᶜ) - Ŝ
    end
    return δx̂ᶜ, PolarFromCartesian(V * (δx̂ᶜ + x̂ᶜ))
end
function NSBPF!(L, U, P, S, x)
    b = S - s(x)
    while norm(b, Inf) > 1e-5
        ### construct d(V) and d(I) from x
        D = dṼ \ dĨ
        z = FEBS(L, U, P, dṼ\b)
        y = z
        for i = 1:5 
            z = FEBS(L, U, P, D*z)
            y += (-1)^i * z 
        end
        x .-= R(Ṽ) \ y
        b = S - s(x)
    end
    return x
end
function DSE!(V, Ĵ, Ĝ, Ĥ, δx̂ᶜ, ŝ₀, s₀, x, Jᶜ, Hᶜ)
    xᶜ  = CartesianFromPolar(x)
    xᵛ   = xᶜ - V * V' * xᶜ
    if norm(xᵛ) > 1e-5
        δx̂ᶜ = push!(δx̂ᶜ, norm(xᵛ))
        xᵛ   ./= norm(xᵛ)
        xᴶ  = Jᶜ * xᵛ
        Ĝ   = [ Ĝ           Ĵ * xᴶ
                xᴶ' * Ĵ     xᴶ' * xᴶ]
        if length(δx̂ᶜ) ≤ n
            Vᵏ  = sort([kron(xᵛ, V), kron(xᵛ, X)])
            Ĥ   = sort([Ĥ           Ĵ' * Hᶜ * Vᵏ
                        xᴶ' * Ĥ     xᴶ' * Hᶜ * Vᵏ])
        end
        ŝ₀  = push!(ŝ₀, xᴶ, s₀)
        Ĵ   = [Ĵ xᴶ]
        V   = [V xᵛ]
    end
    return V, Ĵ, x̂
end