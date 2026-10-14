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
