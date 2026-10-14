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

function _random_edge!(rng, sys, a, b)
    if rand(rng) < 0.5
        return (connect!(sys, Conductor(ThermalConductance(0.01 + 0.09rand(rng))), a, b), :thermal)
    else
        legs = InternalNode(G_a = 0.01 + 0.09rand(rng), G_b = 0.01 + 0.09rand(rng))
        return (connect!(sys, Conductor(VolumeConductance(1e-9 + 9e-9rand(rng)); dissipation = legs), a, b), :piston)
    end
end
