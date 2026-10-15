# Shared helpers for the test suite.

"""
    random_network(rng; nbodies, nreservoirs = 0, extra_edges = 2)

A random connected network of constant-capacity solids joined by thermal
conductors (spanning tree plus `extra_edges` extra conductors, parallel ones
allowed), optionally with fixed-temperature reservoirs attached to random bodies.
Returns `(sys, bodies, reservoirs, conductors, params)`.
"""
function random_network(rng; nbodies, nreservoirs = 0, extra_edges = 2)
    sys = ThermoSystem()
    m = 0.1 .+ 9.9 .* rand(rng, nbodies)          # kg
    c = 100.0 .+ 900.0 .* rand(rng, nbodies)      # J/(kg K)
    T0 = 250.0 .+ 650.0 .* rand(rng, nbodies)     # K
    bodies = [add!(sys, Body(Symbol(:b, i), ConstantCapacitySolid(m[i], c[i]); T = T0[i])) for i in 1:nbodies]
    edges = Tuple{Int,Int}[]
    for i in 2:nbodies
        push!(edges, (rand(rng, 1:i-1), i))
    end
    for _ in 1:extra_edges
        i, j = rand(rng, 1:nbodies), rand(rng, 1:nbodies)
        i == j || push!(edges, (i, j))
    end
    G = 0.1 .+ 9.9 .* rand(rng, length(edges))    # W/K
    conductors = [connect!(sys, Conductor(ThermalConductance(G[k])), bodies[i], bodies[j])
                  for (k, (i, j)) in enumerate(edges)]
    reservoirs = Symbol[]
    for r in 1:nreservoirs
        res = add!(sys, Reservoir(Symbol(:r, r); T = 250.0 + 650.0 * rand(rng)))
        push!(reservoirs, res)
        push!(conductors, connect!(sys, Conductor(ThermalConductance(0.1 + 9.9 * rand(rng))),
                                   res, bodies[rand(rng, 1:nbodies)]))
    end
    return sys, bodies, reservoirs, conductors, (; m, c, T0, edges, G)
end

"""
    slowest_rate(params)

Smallest non-zero eigenvalue of C⁻¹L for an isolated constant-capacity network
(C = diag(m c), L the conductance Laplacian): the slowest relaxation rate [1/s].
"""
function slowest_rate(p)
    n = length(p.m)
    L = zeros(n, n)
    for ((i, j), g) in zip(p.edges, p.G)
        L[i, i] += g; L[j, j] += g; L[i, j] -= g; L[j, i] -= g
    end
    Cmh = Diagonal(1 ./ sqrt.(p.m .* p.c))
    λ = sort(eigvals(Symmetric(Cmh * L * Cmh)))
    return λ[2]
end

# --- slice 2: gases and gas networks ------------------------------------------------------

gas_he() = IdealGas(:He; R = 2077.1, cv = 1.5 * 2077.1)    # monatomic
gas_n2() = IdealGas(:N2; R = 296.8, cv = 2.5 * 296.8)      # diatomic

"Random valid charges for a material (tests of state-function properties)."
random_charges(rng, mat::ConstantCapacitySolid) = charges_at(mat; T = 50 + 1950rand(rng))
random_charges(rng, g::IdealGas) =
    charges_at(g; T = 100 + 1900rand(rng), p = 10.0^(3 + 4rand(rng)), m = 10.0^(-4 + 5rand(rng)))
random_charges(rng, e::ElasticVolume) = charges_at(e; V = e.V0 * (0.5 + rand(rng)))
function random_charges(rng, c::Composite)
    return vcat((random_charges(rng, p) for p in c.parts)...)
end

"""
    random_gas_network(rng; nbodies, nreservoirs = 0, extra_edges = 2)

Ideal-gas bodies (He or N₂) joined by a random spanning tree plus extra edges. Each edge
is a thermal conductor (a diathermal wall) or a damped piston: a volume conductor whose
friction heat goes to an internal node with legs to both gases. Optional reservoirs carry
both T and p. Returns `(sys, bodies, reservoirs, edges)`, where each edge is
`(label, kind)` and kind is `:thermal` or `:piston`.
"""
function random_gas_network(rng; nbodies, nreservoirs = 0, extra_edges = 2)
    sys = ThermoSystem()
    gases = (gas_he(), gas_n2())
    bodies = Symbol[]
    for i in 1:nbodies
        g = gases[rand(rng, 1:2)]
        T = 250 + 500rand(rng); p = 0.5e5 + 2e5rand(rng); V = 1e-3 * (0.5 + rand(rng))
        push!(bodies, add!(sys, Body(Symbol(:g, i), g; T, V, m = p * V / (g.R * T))))
    end
    pairs = [(rand(rng, 1:i-1), i) for i in 2:nbodies]
    for _ in 1:extra_edges
        i, j = rand(rng, 1:nbodies), rand(rng, 1:nbodies)
        i == j || push!(pairs, (i, j))
    end
    edges = Tuple{Symbol,Symbol}[]
    for (i, j) in pairs
        push!(edges, _random_edge!(rng, sys, bodies[i], bodies[j]))
    end
    reservoirs = Symbol[]
    for r in 1:nreservoirs
        res = add!(sys, Reservoir(Symbol(:r, r); T = 250 + 500rand(rng), p = 0.5e5 + 2e5rand(rng)))
        push!(reservoirs, res)
        push!(edges, _random_edge!(rng, sys, res, bodies[rand(rng, 1:nbodies)]))
    end
    return sys, bodies, reservoirs, edges
end

# --- slice 6: mechanical systems ------------------------------------------------------------

"""
    msd_system(; m, k, c, x0, v0 = 0.0, T = 300.0, F = nothing, g = nothing)

A mass on a spring to a fixed frame (`SpringBranch(frame, mass, spring)`: the spring's
extension is the mass's displacement), a viscous damper to the frame with its heat into a
reservoir `:env` (omitted for c = 0), an optional force `Source` `:force` on the mass and an
optional `GravityField` `:gravity`.
"""
function msd_system(; m, k, c, x0, v0 = 0.0, T = 300.0, F = nothing, g = nothing)
    sys = ThermoSystem()
    mass = add!(sys, Body(:mass, InertialMass(m); w = v0))
    spring = add!(sys, Body(:spring, LinearSpring(k); x = x0))
    frame = add!(sys, Reservoir(:frame, Momentum() => 0.0))
    env = add!(sys, Reservoir(:env; T))
    connect!(sys, SpringBranch(), frame, mass, spring; name = :branch)
    c > 0 && connect!(sys, Conductor(LinearConductance(Momentum(), c); dissipation = Into(env)), mass, frame; name = :damper)
    F === nothing || connect!(sys, Source(Momentum() => F), mass; name = :force)
    g === nothing || connect!(sys, GravityField(g), mass; name = :gravity)
    return sys
end

"Free response of a mass–spring–damper released from x0 at rest (08 G1)."
function msd_x(t; ζ, ωn, x0)
    if ζ < 1
        ωd = ωn * sqrt(1 - ζ^2)
        return x0 * exp(-ζ * ωn * t) * (cos(ωd * t) + ζ / sqrt(1 - ζ^2) * sin(ωd * t))
    elseif ζ == 1
        return x0 * exp(-ωn * t) * (1 + ωn * t)
    else
        r = sqrt(ζ^2 - 1)
        s1, s2 = -ωn * (ζ - r), -ωn * (ζ + r)
        return x0 * (s1 * exp(s2 * t) - s2 * exp(s1 * t)) / (s1 - s2)
    end
end

"""
    pendulum_system(; m, ℓ, g, θ0, b = 0.0, T = 300.0)

A pendulum as two accumulators: the rod's inertia mℓ² (`:rod`) and the gravitational potential
mgℓ(1 − cos θ) (`:pot`), joined by a gyrator (dL/dt = −M, dθ/dt = ω). With b > 0, pivot friction
`:pivot` to a frame, its heat into `:env`.
"""
function pendulum_system(; m, ℓ, g, θ0, b = 0.0, T = 300.0)
    sys = ThermoSystem()
    rod = add!(sys, Body(:rod, RotationalInertia(m * ℓ^2); ω = 0.0))
    pot = add!(sys, Body(:pot, PendulumPotential(m, g, ℓ); θ = θ0))
    connect!(sys, Gyrator(-1.0; a = AngularMomentum(), b = Angle()), rod, pot)
    if b > 0
        frame = add!(sys, Reservoir(:frame, AngularMomentum() => 0.0))
        env = add!(sys, Reservoir(:env; T))
        connect!(sys, Conductor(LinearConductance(AngularMomentum(), b); dissipation = Into(env)), rod, frame; name = :pivot)
    end
    return sys
end

"Complete elliptic integral of the first kind K(k), by the arithmetic–geometric mean."
function ellipK(k)
    a, b = 1.0, sqrt(1 - k^2)
    while abs(a - b) > 4eps() * a
        a, b = (a + b) / 2, sqrt(a * b)
    end
    return π / (2a)
end

function _random_edge!(rng, sys, a, b)
    if rand(rng) < 0.5
        return (connect!(sys, Conductor(ThermalConductance(0.01 + 0.09rand(rng))), a, b), :thermal)
    else
        legs = InternalNode(G_a = 0.01 + 0.09rand(rng), G_b = 0.01 + 0.09rand(rng))
        return (connect!(sys, Conductor(VolumeConductance(1e-9 + 9e-9rand(rng)); dissipation = legs), a, b), :piston)
    end
end
