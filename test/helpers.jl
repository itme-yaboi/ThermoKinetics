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

# --- slice 7: species data, mixtures, reactions ----------------------------------------------

const R̄ = 8.314462618                 # J/(mol K)
const N_AV = 6.02214076e23             # 1/mol
const keV = 1.602176634e-16            # J
const YEAR = 365.25 * 86400.0          # s

"Ideal-gas species from molar data at 298.15 K (c_p, h_f, s°), with constant c_v = c_p − R̄."
gasdata(name, atoms; cp, h_f = 0.0, s0 = 0.0) = SpeciesData(name; atoms, cv = cp - R̄, h_f, s0, basis = :molar)
"Condensed species from molar data (constant specific heat c)."
liqdata(name, atoms; c, h_f = 0.0, s0 = 0.0, kw...) = SpeciesData(name; atoms, cv = c, h_f, s0, basis = :molar, kw...)

H2g() = gasdata(:H2, (H = 2,); cp = 28.836, s0 = 130.680)
O2g() = gasdata(:O2, (O = 2,); cp = 29.376, s0 = 205.152)
H2Og() = gasdata(:H2O, (H = 2, O = 1); cp = 33.588, h_f = -241826.0, s0 = 188.835)
N2g() = gasdata(:N2, (N = 2,); cp = 29.124, s0 = 191.609)
I2g() = gasdata(:I2, (I = 2,); cp = 36.888, h_f = 62421.0, s0 = 260.687)
HIg() = gasdata(:HI, (H = 1, I = 1); cp = 29.158, h_f = 26500.0, s0 = 206.594)

"Molar ideal-gas properties from species data (independent of the library's U): u, h, s at (T, p_i)."
struct MolarGas
    cv::Float64; cp::Float64; hf::Float64; s0::Float64
end
MolarGas(sd) = (M = molar_mass(sd); MolarGas(sd.cv * M, sd.cv * M + R̄, sd.h_f * M, sd.s0 * M))
mu_(g::MolarGas, T) = g.hf - R̄ * 298.15 + g.cv * (T - 298.15)
mh_(g::MolarGas, T) = g.hf + g.cp * (T - 298.15)
ms_(g::MolarGas, T, p) = g.s0 + g.cp * log(T / 298.15) - R̄ * log(p / 101325.0)

"Golden-section maximum of a unimodal f on [a, b]."
function golden_max(f, a, b; tol = 1e-15)
    φ = (sqrt(5) - 1) / 2
    c, d = b - φ * (b - a), a + φ * (b - a)
    fc, fd = f(c), f(d)
    while abs(b - a) > tol * max(abs(a), abs(b), 1e-300)
        if fc > fd
            b, d, fd = d, c, fc
            c = b - φ * (b - a); fc = f(c)
        else
            a, c, fc = c, d, fd
            d = a + φ * (b - a); fd = f(d)
        end
    end
    return (a + b) / 2
end

"Moles of each species of a mixture body (NamedTuple by species name)."
function moles(cs, u, body)
    i = findfirst(==(body), body_names(cs))
    mat = cs.bodies[i].material
    E = charges(cs, u, body)
    names = species_names(mat)
    off = length(E) - length(names)
    return NamedTuple{names}(Tuple(E[off+k] / molar_mass(mat, names[k]) for k in eachindex(names)))
end
