# Bodies, reservoirs and conductors (DESIGN §3, §5). Elements hold parameters and laws only;
# they never write state derivatives (that is the assembler's job, assembly.jl).

"""
    Body(name, material, E₀)
    Body(name, material; state...)

A body: a material plus its initial charges. Its charges are its state. The keyword form
passes `state` to `charges_at(material; state...)`, e.g. `Body(:b, gas; T, p, m)`.
"""
struct Body{M<:Material,N}
    name::Symbol
    material::M
    E0::SVector{N,Float64}
end

function Body(name::Symbol, mat::Material, E0::AbstractVector)
    N = ncharges(mat)
    length(E0) == N || throw(ArgumentError("$(typeof(mat)) has $N charges, got $(length(E0))"))
    E = SVector{N,Float64}(E0)
    valid_state(mat, E) || throw(ArgumentError("initial charges $E are not a valid state of $(typeof(mat))"))
    return Body{typeof(mat),N}(name, mat, E)
end
Body(name::Symbol, mat::Material; state...) = Body(name, mat, charges_at(mat; state...))

domains(b::Body) = domains(b.material)

"""
    Reservoir(name; T = nothing, p = nothing, material = nothing)
    Reservoir(name, domain => force, ...)

Fixed forces, infinite capacity (Dirichlet). The keyword form takes temperature and pressure
(stored as the force −p). With a `material` (a fluid with thermal, volumetric and species
charges), the reservoir also fixes the chemical potentials μ(T, p) and the specific entropy
of that fluid, so it can supply or absorb a stream. The pair form takes raw forces,
P = ∂U/∂E, for any domain.
"""
struct Reservoir
    name::Symbol
    domains::Vector{FormOfMotion}
    forces::Vector{Float64}
    material::Union{Nothing,Material}   # fluid held, if the reservoir can exchange mass
    s::Float64                          # specific entropy of that fluid [J/(kg K)] (NaN without one)
    y::Vector{Float64}                  # its mass fractions, in the order of its Species domains
    function Reservoir(name::Symbol, pairs::Pair{<:FormOfMotion,<:Real}...)
        doms = FormOfMotion[first(p) for p in pairs]
        vals = Float64[last(p) for p in pairs]
        return Reservoir(name, doms, vals, nothing, NaN, Float64[])
    end
    function Reservoir(name::Symbol, doms::Vector{FormOfMotion}, vals::Vector{Float64},
                       material::Union{Nothing,Material}, s::Float64, y::Vector{Float64})
        isempty(doms) && throw(ArgumentError("a reservoir needs at least one fixed force"))
        allunique(doms) || throw(ArgumentError("duplicate domains in reservoir :$name"))
        all(isfinite, vals) || throw(ArgumentError("reservoir forces must be finite"))
        for (d, v) in zip(doms, vals)
            d isa Thermal && _check_positive(v, "reservoir temperature")
        end
        return new(name, doms, vals, material, s, y)
    end
end

function Reservoir(name::Symbol; T = nothing, p = nothing, material::Union{Nothing,Material} = nothing)
    if material !== nothing
        (T === nothing || p === nothing) && throw(ArgumentError("a reservoir with a material needs both T and p"))
        _check_positive(T, "reservoir temperature")
        _check_positive(p, "reservoir pressure")
        E = charges_at(material; T, p, m = 1.0)            # one kilogram at (T, p)
        doms = collect(FormOfMotion, domains(material))
        m = sum(E[i] for i in mass_slots(material))
        return Reservoir(name, doms, Vector{Float64}(forces(material, E)), material,
                          E[thermal_slot(material)] / m, Float64[E[i] / m for i in mass_slots(material)])
    end
    pairs = Pair{FormOfMotion,Float64}[]
    T === nothing || (_check_positive(T, "reservoir temperature"); push!(pairs, Thermal() => T))
    p === nothing || (_check_positive(p, "reservoir pressure"); push!(pairs, Volumetric() => -p))
    isempty(pairs) && throw(ArgumentError("give at least one of T or p"))
    return Reservoir(name, pairs...)
end

domains(r::Reservoir) = Tuple(r.domains)

function reservoir_force(r::Reservoir, d::FormOfMotion)
    i = findfirst(==(d), r.domains)
    i === nothing && throw(ArgumentError("reservoir :$(r.name) has no $(d) force"))
    return r.forces[i]
end

has_domain(x::Union{Body,Reservoir}, d::FormOfMotion) = d in domains(x)
has_thermal(x::Union{Body,Reservoir}) = has_domain(x, Thermal())

# --- conductor laws ------------------------------------------------------------------------

abstract type ConductorLaw end

"Whether the law has a thermal row (heat flux J_Q between thermal ports)."
has_thermal_row(::ConductorLaw) = false
"The law's non-thermal rows, as a tuple of domains."
nonthermal_domains(::ConductorLaw) = ()

"""
    ThermalConductance(G)

Energy-form thermal conductance G ≥ 0 [W/K]: Q̇ = G (T_a − T_b). Implemented through the exact
discrete kernel A (07 §3.1) with M_QQ = G T_a T_b and X_Q = 1/T_b − 1/T_a, so the entropy
production Q̇ (1/T_b − 1/T_a) is ≥ 0 at any ΔT. The sign of G is a type invariant.
"""
struct ThermalConductance <: ConductorLaw
    G::Float64
    function ThermalConductance(G::Real)
        (G >= 0 && isfinite(G)) || throw(ArgumentError("conductance must be ≥ 0 and finite, got $G"))
        return new(G)
    end
end

has_thermal_row(::ThermalConductance) = true

"Kernel A coefficient M_QQ for the thermal row at end temperatures `Ta`, `Tb`."
conductance_QQ(law::ThermalConductance, Ta, Tb) = law.G * Ta * Tb

"""
    VolumeConductance(L)

A damped piston or flexible membrane between two volumetric ports (07 §5): volume flux
J_V = L (P_a − P_b) with P = −p, L ≥ 0 [m³/(Pa s)]. Volume moves toward the higher pressure,
which is correct for a moving boundary. Pressure-driven *fluid* flow is mass transport,
not this. The released power D = L ΔP² ≥ 0 is friction heat, placed by the conductor's
dissipation policy.
"""
struct VolumeConductance <: ConductorLaw
    L::Float64
    function VolumeConductance(L::Real)
        (L >= 0 && isfinite(L)) || throw(ArgumentError("volume conductance must be ≥ 0 and finite, got $L"))
        return new(L)
    end
end

nonthermal_domains(::VolumeConductance) = (Volumetric(),)
fluxes(law::VolumeConductance, ΔP) = SVector(law.L * ΔP[1])

# --- dissipation policies (DESIGN §5.1) ------------------------------------------------------

abstract type DissipationPolicy end

"""
    Split(w)

Kernel A: a fraction `w` of the released power goes to port b as heat, 1 − w to port a.
"""
struct Split <: DissipationPolicy
    w::Float64
    function Split(w::Real)
        (0 <= w <= 1) || throw(ArgumentError("Split fraction must be in [0, 1], got $w"))
        return new(w)
    end
end

"""
    Into(node)

Kernel A: all released power goes as heat to a third thermal node (a wall or a reservoir).
"""
struct Into <: DissipationPolicy
    node::Symbol
end

"""
    InternalNode(; G_a = 0, G_b = 0, G_h = 0, h = nothing)

Kernel B, the default for non-thermal dissipation: a zero-capacity thermal node receives the
released power and conducts it through legs of conductance G [W/K] to port a, port b and an
optional third node `h`. Its temperature is eliminated analytically,
T_m = (Σ G_i T_i + D)/Σ G_i, so the model stays an ODE.
"""
struct InternalNode <: DissipationPolicy
    G_a::Float64
    G_b::Float64
    G_h::Float64
    h::Union{Symbol,Nothing}
    function InternalNode(; G_a::Real = 0.0, G_b::Real = 0.0, G_h::Real = 0.0, h::Union{Symbol,Nothing} = nothing)
        for (n, G) in (("G_a", G_a), ("G_b", G_b), ("G_h", G_h))
            (G >= 0 && isfinite(G)) || throw(ArgumentError("$n must be ≥ 0 and finite, got $G"))
        end
        G_a + G_b + G_h > 0 || throw(ArgumentError("an internal node needs at least one leg with G > 0"))
        (G_h > 0) == (h !== nothing) || throw(ArgumentError("give G_h > 0 together with the third node h, or neither"))
        return new(G_a, G_b, G_h, h)
    end
end

third_node(::Union{Nothing,Split}) = nothing
third_node(p::Into) = p.node
third_node(p::InternalNode) = p.h

"""
    Conductor(law; dissipation = nothing)

A conductor element, joined to two nodes with `connect!`. A law with non-thermal rows
releases power D ≥ 0, so it must say where that heat goes: `InternalNode(G_a = …, G_b = …)`
(the default kind, kernel B), `Split(w)` or `Into(node)` (kernel A). A thermal-only law takes
no policy.
"""
struct Conductor{L<:ConductorLaw,P<:Union{Nothing,DissipationPolicy}}
    law::L
    dissipation::P
    function Conductor(law::ConductorLaw; dissipation::Union{Nothing,DissipationPolicy} = nothing)
        rows = nonthermal_domains(law)
        if isempty(rows)
            dissipation === nothing ||
                throw(ArgumentError("$(typeof(law)) has no non-thermal rows; it takes no dissipation policy"))
        else
            has_thermal_row(law) &&
                throw(ArgumentError("coupled conductors (thermal + non-thermal rows) arrive in slice 4"))
            dissipation === nothing &&
                throw(ArgumentError("$(typeof(law)) dissipates; give dissipation = InternalNode(G_a = …, G_b = …) " *
                                    "(default kind), Split(w) or Into(node)"))
        end
        return new{typeof(law),typeof(dissipation)}(law, dissipation)
    end
end

domains(c::Conductor) = has_thermal_row(c.law) ? (Thermal(), nonthermal_domains(c.law)...) : nonthermal_domains(c.law)

# --- streams, pumps, flow boundaries (DESIGN §7) ---------------------------------------------

abstract type FlowDriver end

"""
    HydraulicConductance(κ)

Flow driver for a `Stream`: ṁ = κ (p_a − p_b), κ ≥ 0 [kg/(s Pa)]. The sign of ṁ follows the
pressure difference by type, so a stream between bodies always runs down the pressure
gradient. A driver only sets ṁ: it deposits no friction heat (the stream's own generation
already contains the throttling loss, 07 §4).
"""
struct HydraulicConductance <: FlowDriver
    κ::Float64
    function HydraulicConductance(κ::Real)
        (κ >= 0 && isfinite(κ)) || throw(ArgumentError("hydraulic conductance must be ≥ 0 and finite, got $κ"))
        return new(κ)
    end
end

mass_flow(d::HydraulicConductance, pa, pb) = d.κ * (pa - pb)

"""
    Stream(driver)

Advection of a fluid between two nodes (07 §4): a bundle of charges, namely each species'
mass and the entropy ṁ s, evaluated at the upstream state (upwinding). Volume is not
advected into a fixed control volume. The downstream node receives the mass at its own μ,
and the rest of the carried energy ṁ h_u as entropy at its own T. The generation
ṁ[s_u(T_u − T_d) + Σ y(μ_u − μ_d)]/T_d is booked on the stream.
"""
struct Stream{D<:FlowDriver}
    driver::D
end

"""
    IdealPump(mdot)

Work-supplied flow from port a to port b at a prescribed ṁ = mdot(t) ≥ 0 [kg/s], whatever
the pressures (D-13). The fluid is raised isentropically from p_a to p_b: the delivered state
is (s_a, p_b), and the shaft work ṁ (h(s_a, p_b) − h_a) is booked as external work. Negative
work, with flow down the gradient, is an ideal turbine.
"""
struct IdealPump{F}
    mdot::F
end

"""
    FlowBoundary(material; T, p, mdot, supply)

External inflow of `material` at (T, p) into a node, at ṁ = mdot(t) ≥ 0 [kg/s]. `supply` is
required:
- `:passive`: the fluid arrives at (T, p). The node pressure must not exceed p (an error
  otherwise, naming the boundary).
- `:pumped`: an ideal pump raises it isentropically to the node pressure, with the work booked.

The advected entropy is booked as boundary flow; the inlet-face generation is booked on the
boundary.
"""
struct FlowBoundary{M<:Material,F}
    material::M
    T::Float64
    p::Float64
    mdot::F
    supply::Symbol
    function FlowBoundary(material::Material; T::Real, p::Real, mdot, supply::Union{Symbol,Nothing} = nothing)
        supply === nothing &&
            throw(ArgumentError("label the flow boundary: supply = :passive (fluid arrives at its own p) or :pumped"))
        supply in (:passive, :pumped) || throw(ArgumentError("supply must be :passive or :pumped, got :$supply"))
        _check_positive(T, "boundary temperature")
        _check_positive(p, "boundary pressure")
        return new{typeof(material),typeof(mdot)}(material, T, p, mdot, supply)
    end
end

"""
    PassThroughOutlet(; mdot)
    PassThroughOutlet(; κ, p_out)

External outflow (07 §5): carries the node's own bundle out of the model with no
generation. Either a prescribed ṁ = mdot(t) ≥ 0, or ṁ = κ max(p − p_out, 0) to a back
pressure (no backflow: an outlet cannot supply a state).
"""
struct PassThroughOutlet{F}
    mdot::F                 # a function of t, or `nothing`
    κ::Float64
    p_out::Float64
    function PassThroughOutlet(; mdot = nothing, κ = nothing, p_out = nothing)
        if mdot !== nothing
            (κ === nothing && p_out === nothing) || throw(ArgumentError("give either mdot, or κ and p_out"))
            return new{typeof(mdot)}(mdot, NaN, NaN)
        end
        (κ === nothing || p_out === nothing) && throw(ArgumentError("give either mdot, or κ and p_out"))
        (κ >= 0 && isfinite(κ)) || throw(ArgumentError("outlet conductance must be ≥ 0 and finite, got $κ"))
        _check_positive(p_out, "back pressure")
        return new{Nothing}(nothing, Float64(κ), Float64(p_out))
    end
end

outlet_flow(o::PassThroughOutlet{Nothing}, p, t) = o.κ * max(p - o.p_out, zero(p))
outlet_flow(o::PassThroughOutlet, p, t) = _nonneg_flow(o.mdot(t), "outlet")

function _nonneg_flow(m, what)
    m >= 0 || throw(ArgumentError("prescribed $what mass flow must be ≥ 0, got $m"))
    return m
end
