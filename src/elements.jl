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
P = ∂U/∂E, for any domain. A force may be a function of time, `domain => t -> P(t)`: a
prescribed boundary history, such as a drive velocity. Pass its discontinuities to
`simulate` as `tstops`.
"""
struct Reservoir
    name::Symbol
    domains::Vector{FormOfMotion}
    forces::Vector{Float64}             # values at t = 0
    fns::Vector{Any}                    # `nothing`, or a function t -> force, per domain
    material::Union{Nothing,Material}   # fluid held, if the reservoir can exchange mass
    s::Float64                          # specific entropy of that fluid [J/(kg K)] (NaN without one)
    y::Vector{Float64}                  # its mass fractions, in the order of its Species domains
    function Reservoir(name::Symbol, pairs::Pair{<:FormOfMotion}...)
        doms = FormOfMotion[first(p) for p in pairs]
        fns = Any[last(p) isa Real ? nothing : last(p) for p in pairs]
        vals = Float64[last(p) isa Real ? last(p) : last(p)(0.0) for p in pairs]
        return Reservoir(name, doms, vals, nothing, NaN, Float64[]; fns)
    end
    function Reservoir(name::Symbol, doms::Vector{FormOfMotion}, vals::Vector{Float64},
                       material::Union{Nothing,Material}, s::Float64, y::Vector{Float64};
                       fns::Vector{Any} = Any[nothing for _ in doms])
        isempty(doms) && throw(ArgumentError("a reservoir needs at least one fixed force"))
        allunique(doms) || throw(ArgumentError("duplicate domains in reservoir :$name"))
        all(isfinite, vals) || throw(ArgumentError("reservoir forces must be finite"))
        for (d, v) in zip(doms, vals)
            d isa Thermal && _check_positive(v, "reservoir temperature")
        end
        return new(name, doms, vals, fns, material, s, y)
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

function reservoir_force(r::Reservoir, d::FormOfMotion, t::Real = 0.0)
    i = findfirst(==(d), r.domains)
    i === nothing && throw(ArgumentError("reservoir :$(r.name) has no $(d) force"))
    return r.fns[i] === nothing ? r.forces[i] : Float64(r.fns[i](t))
end

"""
    Junction(name, domains...)

A zero-capacity node: it holds no charge and no energy. Its forces are algebraic unknowns
fixed by "net charge flux in = 0" for each domain (a Kirchhoff current law, a massless
mechanical node): an index-1 DAE row in the mass-matrix form (DESIGN §6.2). Thermal and
species junctions are not allowed: a zero-capacity node cannot hold entropy, and mixing needs
a thermal DOF.
"""
struct Junction
    name::Symbol
    domains::Vector{FormOfMotion}
    function Junction(name::Symbol, doms::FormOfMotion...)
        isempty(doms) && throw(ArgumentError("a junction needs at least one form of motion"))
        allunique(doms) || throw(ArgumentError("duplicate domains in junction :$name"))
        any(d -> d isa Thermal || d isa Species, doms) &&
            throw(ArgumentError("junctions are non-thermal and hold no mass: no Thermal or Species"))
        return new(name, collect(FormOfMotion, doms))
    end
end

domains(j::Junction) = Tuple(j.domains)

has_domain(x::Union{Body,Reservoir,Junction}, d::FormOfMotion) = d in domains(x)
has_thermal(x::Union{Body,Reservoir,Junction}) = has_domain(x, Thermal())

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
    LinearConductance(domain, L)

One non-thermal row: J = L (P_a − P_b), L ≥ 0, for any form of motion. The released power
D = L ΔP² ≥ 0 is placed by the conductor's dissipation policy. Use the named constructors
`VolumeConductance`, `Resistor` and `MassConductance`.
"""
struct LinearConductance{D<:FormOfMotion} <: ConductorLaw
    domain::D
    L::Float64
    function LinearConductance(domain::FormOfMotion, L::Real)
        domain isa Thermal && throw(ArgumentError("use ThermalConductance for the thermal row"))
        (L >= 0 && isfinite(L)) || throw(ArgumentError("conductance must be ≥ 0 and finite, got $L"))
        return new{typeof(domain)}(domain, L)
    end
end

nonthermal_domains(law::LinearConductance) = (law.domain,)
fluxes(law::LinearConductance, ΔP) = SVector(law.L * ΔP[1])

"""
    VolumeConductance(L)

A damped piston or flexible membrane between two volumetric ports (07 §5): volume flux
J_V = L (P_a − P_b) with P = −p, L ≥ 0 [m³/(Pa s)]. Volume moves toward the higher pressure,
which is correct for a moving boundary. Pressure-driven *fluid* flow is mass transport,
not this.
"""
VolumeConductance(L::Real) = LinearConductance(Volumetric(), L)

"""
    Resistor(R)

Ohmic conductor, R > 0 [Ω]: J_e = (φ_a − φ_b)/R. Joule heat D = ΔV²/R is placed by the
dissipation policy (`InternalNode` for a resistor with its own heat path, `Into(node)` to
dump it into a wall or reservoir).
"""
function Resistor(R::Real)
    (R > 0 && isfinite(R)) || throw(ArgumentError("resistance must be positive and finite, got $R"))
    return LinearConductance(Electrical(), 1 / R)
end

"""
    MassConductance(species, L)

Diffusion of one species driven by its chemical potential: J_m = L (μ_a − μ_b), L ≥ 0
[kg² /(J s)]. Using μ, not concentration, is the contact rule (non-negotiable 5): at
equilibrium μ is equal across the barrier, while p and concentration need not be.
"""
MassConductance(species::Symbol, L::Real) = LinearConductance(Species(species), L)

# --- coupled laws (a thermal row plus non-thermal rows; kernel A, 07 §3.1) -------------------
#
# The forces are x = (X_Q, Y_1…Y_N) with X_Q = 1/T_b − 1/T_a and Y_k = ΔP_k c,
# c = w_a/T_a + w_b/T_b. A coupled law supplies a factor B and an antisymmetric W at the
# current state: J = (BBᵀ + W) x. Entropy production is ‖Bᵀx‖² ≥ 0 by construction.

"""
    coupled_factor(law, Ta, Tb, Pa, Pb, c, w) -> (B, W)

Factor B (M_sym = BBᵀ) and antisymmetric part W of a coupled law at end temperatures
`Ta`, `Tb`, end forces `Pa`, `Pb` of its non-thermal rows, weight `c` and Split fraction `w`.
"""
function coupled_factor end

"The full coupled matrix BBᵀ + W at a state, for inspection and tests (c from the Split fraction w)."
function coupled_matrix(law::ConductorLaw, Ta, Tb, Pa, Pb, w)
    c = (1 - w) / Ta + w / Tb
    B, W = coupled_factor(law, Ta, Tb, Pa, Pb, c, w)
    return B * B' + W
end

"Two-row factor B = [√K′ √g·b; 0 √g]: J_k = g (Y + b X_Q), J_Q = g b Y + (g b² + K′) X_Q."
@inline function _two_row_factor(Kp, g, b, Kzero::Bool)
    sK = Kzero ? zero(g) : sqrt(Kp)
    sg = sqrt(g)
    return SMatrix{2,2}(sK, zero(sg), sg * b, sg), zero(SMatrix{2,2,typeof(sg)})
end

"""
    ThermoelectricConductance(; R, K, α)

A thermoelectric element (Veinik's thermoelectric couple leg, §44–47): electrical resistance
R > 0 [Ω], open-circuit thermal conductance K ≥ 0 [W/K], Seebeck coefficient α [V/K].
Calibrated so that, exactly at any finite ΔT:
- J_e = (Δφ + α ΔT)/R, so the open-circuit voltage is −α ΔT (Seebeck; (278));
- J_Q = K ΔT + α T_a T_b c J_e, the Peltier heat Π J_e with Π = α T̄ (Kelvin's second relation);
- the coupled matrix is PSD with margin exactly K T_a T_b.
"""
struct ThermoelectricConductance <: ConductorLaw
    R::Float64
    K::Float64
    α::Float64
    function ThermoelectricConductance(; R::Real, K::Real, α::Real)
        (R > 0 && isfinite(R)) || throw(ArgumentError("resistance must be positive and finite, got $R"))
        (K >= 0 && isfinite(K)) || throw(ArgumentError("thermal conductance must be ≥ 0 and finite, got $K"))
        isfinite(α) || throw(ArgumentError("Seebeck coefficient must be finite"))
        return new(R, K, α)
    end
end

has_thermal_row(::ThermoelectricConductance) = true
nonthermal_domains(::ThermoelectricConductance) = (Electrical(),)

function coupled_factor(law::ThermoelectricConductance, Ta, Tb, Pa, Pb, c, w)
    g = 1 / (law.R * c)
    b = law.α * Ta * Tb * c
    return _two_row_factor(law.K * Ta * Tb, g, b, law.K == 0)
end

"""
    thermoelectric_from_veinik(L11, L12, L22; T)

Convert a leg's coefficients in Veinik's (277) variables at temperature T into a
`ThermoelectricConductance`:
- J_E = L11 δφ + L12 ΔT/T;
- J_Q = L21 δφ + L22 ΔT/T.

R = 1/L11, α = L12/(L11 T), K = L22/T − α² T L11. Requires L11 L22 ≥ L12² (that is, K ≥ 0).
"""
function thermoelectric_from_veinik(L11::Real, L12::Real, L22::Real; T::Real)
    _check_positive(L11, "L11")
    _check_positive(T, "temperature")
    α = L12 / (L11 * T)
    K = L22 / T - α^2 * T * L11
    K >= -1e-12 * abs(L22 / T) ||
        throw(ArgumentError("L11 L22 < L12²: the leg's Onsager matrix is not positive semidefinite"))
    return ThermoelectricConductance(R = 1 / L11, K = max(K, 0.0), α = α)
end

"""
    KnudsenBarrier(gas; g, K = 0.0)

A molecular-flow (Knudsen) barrier between two volumes of `gas` (an `IdealGas`), with mass and
heat rows (Veinik §49). Mass flux J_m = g F, with F = Δ(μ/T) + (c_p − R/2) ln(T_a/T_b) and
g > 0 [kg² K/(J s)]. Hence J_m = 0 exactly when p_a/√T_a = p_b/√T_b (thermal transpiration,
(300)), at any finite ΔT. K ≥ 0 [W/K] is the barrier's own heat conduction. Near equilibrium
it coincides with Veinik's barrier (291)–(299) transformed to the mass basis
(`barrier_mass_basis`, with ξ = g v²/T).
"""
struct KnudsenBarrier <: ConductorLaw
    gas::IdealGas
    g::Float64
    K::Float64
    function KnudsenBarrier(gas::IdealGas; g::Real, K::Real = 0.0)
        (g > 0 && isfinite(g)) || throw(ArgumentError("g must be positive and finite, got $g"))
        (K >= 0 && isfinite(K)) || throw(ArgumentError("thermal conductance must be ≥ 0 and finite, got $K"))
        return new(gas, g, K)
    end
end

has_thermal_row(::KnudsenBarrier) = true
nonthermal_domains(law::KnudsenBarrier) = (Species(law.gas.species),)

"ln(r)/(r − 1), smooth through r = 1."
@inline function _logratio(r)
    d = r - 1
    return abs(d) < 1e-4 ? 1 - d / 2 + d^2 / 3 - d^3 / 4 : log(r) / d
end

function coupled_factor(law::KnudsenBarrier, Ta, Tb, Pa, Pb, c, w)
    gas = law.gas
    cp = gas.cv + gas.R
    μ̃ = w * Pa[1] + (1 - w) * Pb[1]                      # μ_a/T_a − μ_b/T_b = Δμ c − μ̃ X_Q
    λ = Ta * _logratio(Ta / Tb)                          # ln(T_a/T_b) = λ X_Q
    b = -μ̃ + (cp - gas.R / 2) * λ
    return _two_row_factor(law.K * Ta * Tb, law.g + zero(b), b, law.K == 0)
end

"""
    OnsagerConductance(rows, M; W = nothing, magnetic = false)

A general linear conductor over `rows` (a tuple of domains). Validated at construction: M must
be symmetric and positive semidefinite; an antisymmetric W (Hall-type effects) is accepted
only with `magnetic = true` (Onsager–Casimir, non-negotiable 4).
- If `rows` starts with `Thermal()`: a coupled conductor (kernel A). M acts on the exact
  discrete forces x = (X_Q, ΔP c) and is constant; build it from data with
  `barrier_mass_basis` or similar.
- Otherwise: non-thermal rows only, J = (M + W) ΔP, with heat placed by the dissipation
  policy.
"""
struct OnsagerConductance{R<:Tuple,N,NN} <: ConductorLaw
    rows::R
    M::SMatrix{N,N,Float64,NN}
    B::SMatrix{N,N,Float64,NN}
    W::SMatrix{N,N,Float64,NN}
    function OnsagerConductance(rows::Tuple, M::AbstractMatrix; W = nothing, magnetic::Bool = false)
        N = length(rows)
        size(M) == (N, N) || throw(ArgumentError("M must be $N × $N for $N rows"))
        all(d -> d isa FormOfMotion, rows) && allunique(rows) || throw(ArgumentError("rows must be distinct domains"))
        nth = count(d -> d isa Thermal, rows)
        (nth == 0 || (nth == 1 && rows[1] isa Thermal)) ||
            throw(ArgumentError("the thermal row, if any, must come first"))
        scale = max(maximum(abs, M), floatmin())
        maximum(abs, M - M') <= 1e-12 * scale ||
            throw(ArgumentError("M must be symmetric (Onsager); an antisymmetric part goes in W, with magnetic = true"))
        F = eigen(Symmetric(Matrix{Float64}(M)))
        minimum(F.values) >= -1e-12 * scale ||
            throw(ArgumentError("M must be positive semidefinite (second law); eigenvalues $(F.values)"))
        B = F.vectors * Diagonal(sqrt.(max.(F.values, 0.0)))
        Wm = zeros(N, N)
        if W !== nothing
            magnetic || throw(ArgumentError("an antisymmetric part W is only allowed in a magnetic field: magnetic = true"))
            size(W) == (N, N) && W == -transpose(W) || throw(ArgumentError("W must be exactly antisymmetric"))
            Wm = Matrix{Float64}(W)
        end
        return new{typeof(rows),N,N * N}(rows, SMatrix{N,N}(Symmetric(Matrix{Float64}(M))), SMatrix{N,N}(B), SMatrix{N,N}(Wm))
    end
end

has_thermal_row(law::OnsagerConductance) = law.rows[1] isa Thermal
nonthermal_domains(law::OnsagerConductance) = has_thermal_row(law) ? Base.tail(law.rows) : law.rows
fluxes(law::OnsagerConductance, ΔP) = (law.M + law.W) * ΔP
coupled_factor(law::OnsagerConductance, Ta, Tb, Pa, Pb, c, w) = (law.B, law.W)

"""
    barrier_mass_basis(LV; T, v, s) -> M

Transform Veinik's barrier coefficients LV, in the (J_v, J_s) basis against (Δp, ΔT) with
Tσ = J·X ((291), lumped), into the library's coupled (heat, mass) basis at the reference state
(T, v, s) of the fluid. The result is the congruence M = T·A·LV·Aᵀ with
A = [1/v 0; T s/v T] in (mass, heat) order, returned in (heat, mass) order. Symmetry and
positive semidefiniteness are preserved, and so is entropy production: σ = J_V·X_V/T = Jᵀx.
"""
function barrier_mass_basis(LV::AbstractMatrix; T::Real, v::Real, s::Real)
    size(LV) == (2, 2) || throw(ArgumentError("LV must be 2 × 2"))
    A = [1/v 0.0; T*s/v T]
    Mmq = T * A * LV * A'
    return [Mmq[2, 2] Mmq[2, 1]; Mmq[1, 2] Mmq[1, 1]]
end

"""
    electrokinetic_mass_basis(LV; v) -> M

Transform Veinik's electrokinetic coefficients LV, in the (J_e, J_v) basis against (Δφ, Δp)
((316)–(317), isothermal), into the (charge, mass) basis against (Δφ, Δμ), with Δμ = v Δp at
fixed T: M = A·LV·Aᵀ with A = diag(1, 1/v).
"""
function electrokinetic_mass_basis(LV::AbstractMatrix; v::Real)
    size(LV) == (2, 2) || throw(ArgumentError("LV must be 2 × 2"))
    A = [1.0 0.0; 0.0 1/v]
    return A * LV * A'
end

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

A conductor element, joined to two nodes with `connect!`.
- A thermal-only law takes no policy.
- A law with non-thermal rows only releases power D ≥ 0, so it must say where that heat
  goes: `InternalNode(G_a = …, G_b = …)` (the default kind, kernel B), `Split(w)` or
  `Into(node)` (kernel A).
- A coupled law (a thermal row plus non-thermal rows) puts its released power into its own
  ends: `Split(w)` only, defaulting to `Split(0.5)`. With coupling, D can be negative
  (thermoelectric generation), so it may not be routed to a third node.
"""
struct Conductor{L<:ConductorLaw,P<:Union{Nothing,DissipationPolicy}}
    law::L
    dissipation::P
    function Conductor(law::ConductorLaw; dissipation::Union{Nothing,DissipationPolicy} = nothing)
        rows = nonthermal_domains(law)
        if isempty(rows)
            dissipation === nothing ||
                throw(ArgumentError("$(typeof(law)) has no non-thermal rows; it takes no dissipation policy"))
        elseif has_thermal_row(law)
            dissipation === nothing && (dissipation = Split(0.5))
            dissipation isa Split ||
                throw(ArgumentError("a coupled conductor puts its released power into its own ends: use Split(w)"))
        else
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


# --- converters (DESIGN §6) -------------------------------------------------------------------

"""
    Converter(W, domains)

An ideal multi-port converter. Port k, on `domains[k]`, receives the charge flux
F_k = Σ_j W_kj P_j, where P_j are the port forces. W must be exactly antisymmetric, so
Σ P_k F_k = Pᵀ W P ≡ 0: power is conserved and no entropy is produced (the reversible, Poisson
part of GENERIC). It is the one connection that does not conserve each charge
(non-negotiable 1). The thermal form of motion is never converted. Join with
`connect!(sys, conv, node_1, …, node_n)`.
"""
struct Converter{N,NN}
    W::SMatrix{N,N,Float64,NN}
    domains::NTuple{N,FormOfMotion}
    function Converter(W::AbstractMatrix, doms::Tuple)
        N = length(doms)
        N >= 2 || throw(ArgumentError("a converter needs at least two ports"))
        size(W) == (N, N) || throw(ArgumentError("W must be $N × $N for $N ports"))
        all(d -> d isa FormOfMotion, doms) || throw(ArgumentError("ports must be forms of motion"))
        any(d -> d isa Thermal, doms) && throw(ArgumentError("thermal charge (entropy) is never converted"))
        W == -transpose(W) || throw(ArgumentError("W must be exactly antisymmetric (power conservation)"))
        return new{N,N * N}(SMatrix{N,N,Float64}(W), Tuple(doms))
    end
end

"""
    Gyrator(r; a, b)

Two-port converter: J_a = r P_b and J_b = −r P_a. Examples:
- a piston of area A between a gas volume and a mass: `Gyrator(A; a = Volumetric(), b = Momentum())`;
- the atmosphere on the other face: `Gyrator(-A; …)`.
"""
function Gyrator(r::Real; a::FormOfMotion, b::FormOfMotion)
    (isfinite(r) && r != 0) || throw(ArgumentError("gyrator ratio must be finite and non-zero"))
    return Converter([0.0 r; -r 0.0], (a, b))
end

"""
    InductorBranch()

The 3-port converter that puts an inductor body (flux linkage Λ, force i) between two electrical
nodes a and b: dq_a = −i, dq_b = +i, dΛ = φ_a − φ_b. Connect it as
`connect!(sys, InductorBranch(), a, b, inductor)`.
"""
InductorBranch() = Converter([0.0 0.0 -1.0; 0.0 0.0 1.0; 1.0 -1.0 0.0], (Electrical(), Electrical(), FluxLinkage()))

"""
    Transformer(r; a, b, τ = 1e-4)

Ideal transformer between ports a and b: forces P_a = r P_b, flows F_a = λ and F_b = −r λ, so
the power λ (P_a − r P_b) is zero on the constraint. Between two storing bodies this constraint
is index 2. The library uses the index-reduced form: λ is an algebraic unknown whose row is the
differentiated constraint with Baumgarte stabilization, d/dt(P_a − r P_b) + (P_a − r P_b)/τ = 0
(index 1, D-29). The initial state must satisfy the constraint. Ends are bodies or
constant-force reservoirs.
"""
struct Transformer
    r::Float64
    a::FormOfMotion
    b::FormOfMotion
    τ::Float64
    function Transformer(r::Real; a::FormOfMotion, b::FormOfMotion, τ::Real = 1e-4)
        (isfinite(r) && r != 0) || throw(ArgumentError("transformer ratio must be finite and non-zero"))
        (a isa Thermal || b isa Thermal) && throw(ArgumentError("thermal charge (entropy) is never converted"))
        _check_positive(τ, "stabilization time τ")
        return new(Float64(r), a, b, Float64(τ))
    end
end

"""
    Switched(element; on)

Gate an element by a schedule `on(t) -> Bool` (D-14): while off, the element moves nothing and
produces nothing. Pass the switching times to `simulate` as `tstops`. Transformers cannot be
switched (their algebraic row would be undefined).
"""
struct Switched{E,F}
    element::E
    on::F
    function Switched(element; on)
        element isa Transformer && throw(ArgumentError("a transformer cannot be switched"))
        return new{typeof(element),typeof(on)}(element, on)
    end
end
