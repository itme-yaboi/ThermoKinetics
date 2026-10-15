# Internal transformations and sources (DESIGN §8, §8.1; D-11, D-39): reactions inside one
# body (affinity-driven, dU = 0, production A·ξ̇/T ≥ 0), nuclear processes (prescribed rates,
# nucleons and charge conserved, dU = E_dep ξ̇) and heat sources. Elements hold laws only; the
# assembler does the bookkeeping.

# --- reaction laws --------------------------------------------------------------------------------

abstract type ReactionLaw end

"""
    AffinityLaw(L)

Linear affinity law for R reactions in one body: ξ̇ = L·(A/T), with L ≥ 0 [mol² K/(J s)] (a
number for one reaction) or an R × R symmetric positive semidefinite matrix (Onsager coupling
between reactions). The production (A/T)ᵀ L (A/T) is computed as ‖Bᵀ(A/T)‖² with L = B Bᵀ, so
it is ≥ 0 in floating point.
"""
struct AffinityLaw{R,RR} <: ReactionLaw
    L::SMatrix{R,R,Float64,RR}
    B::SMatrix{R,R,Float64,RR}
end

function AffinityLaw(L::Real)
    (L >= 0 && isfinite(L)) || throw(ArgumentError("the affinity coefficient must be ≥ 0 and finite, got $L"))
    return AffinityLaw{1,1}(SMatrix{1,1}(Float64(L)), SMatrix{1,1}(sqrt(Float64(L))))
end

function AffinityLaw(L::AbstractMatrix)
    R = size(L, 1)
    size(L) == (R, R) || throw(ArgumentError("L must be square"))
    all(isfinite, L) || throw(ArgumentError("L must be finite"))
    scale = max(maximum(abs, L), floatmin())
    maximum(abs, L - L') <= 1e-12 * scale ||
        throw(ArgumentError("L must be symmetric (Onsager reciprocity; no magnetic field acts on a reaction)"))
    F = eigen(Symmetric(Matrix{Float64}(L)))
    minimum(F.values) >= -1e-12 * scale ||
        throw(ArgumentError("L must be positive semidefinite (second law); eigenvalues $(F.values)"))
    B = F.vectors * Diagonal(sqrt.(max.(F.values, 0.0)))
    return AffinityLaw{R,R * R}(SMatrix{R,R}(Symmetric(Matrix{Float64}(L))), SMatrix{R,R}(B))
end

"""
    MassAction(k0; b = 0, Ea = 0, orders = nothing)

Mass-action law for one reaction: ξ̇ = V k(T) Π_r c_r^(order_r) · (1 − e^(−A/(R̄T))), with
c = n/V [mol/m³] and k(T) = k0 (T/298.15)^b exp(−Ea/(R̄T)). `orders` (a NamedTuple by
species) defaults to the reactants' stoichiometric coefficients. The reverse rate is implied
by the thermodynamics (detailed balance), so no equilibrium constant is entered. The factor
(1 − e^(−A/R̄T)) is evaluated with `expm1` and has the sign of A exactly, so A·ξ̇ ≥ 0 in
floating point. Needs a body with a volume (`Volumetric`).
"""
struct MassAction{O} <: ReactionLaw
    k0::Float64
    b::Float64
    Ea::Float64
    orders::O
    function MassAction(k0::Real; b::Real = 0.0, Ea::Real = 0.0, orders = nothing)
        (k0 >= 0 && isfinite(k0)) || throw(ArgumentError("the rate constant must be ≥ 0 and finite, got $k0"))
        (isfinite(b) && isfinite(Ea)) || throw(ArgumentError("b and Ea must be finite"))
        orders === nothing || orders isa NamedTuple || throw(ArgumentError("orders must be a NamedTuple by species"))
        orders === nothing || all(o -> o >= 0 && isfinite(o), values(orders)) ||
            throw(ArgumentError("reaction orders must be ≥ 0"))
        return new{typeof(orders)}(Float64(k0), Float64(b), Float64(Ea), orders)
    end
end

_rate_constant(law::MassAction, T) = law.k0 * (T / 298.15)^law.b * exp(-law.Ea / (R_MOLAR * T))

"""
    extent_rate(law::MassAction, A, T, V, n; orders = ones) -> ξ̇

The extent rate [mol/s] at affinity A [J/mol], temperature T, volume V and the amounts `n`
[mol] of the species in the forward rate (with their `orders`). Evaluated in log form, so a
large |A| neither overflows nor loses the sign.
"""
extent_rate(law::MassAction, A, T, V, n::Tuple; orders = map(_ -> 1.0, n)) =
    _mass_action_rate(law, A, T, V, map(log, n), orders)

"Mass action from the log amounts (the activity amounts in a body, D-41)."
function _mass_action_rate(law::MassAction, A, T, V, logn::Tuple, orders)
    k = _rate_constant(law, T)
    lnf = sum(ntuple(j -> orders[j] == 0 ? zero(A) : orders[j] * (logn[j] - log(V)), length(logn)); init = zero(A))
    x = -A / (R_MOLAR * T)
    return x <= 0 ? V * k * exp(lnf) * (-expm1(x)) : -V * k * exp(lnf + x) * (-expm1(-x))
end
extent_rate(law::AffinityLaw{1}, A, T, V, n::Tuple) = law.L[1, 1] * A / T

# --- reactions ------------------------------------------------------------------------------------

"""
    Reaction(ν::Pair...; law)
    Reaction([ν₁, ν₂, …]; law)

Chemical reactions inside one body (DESIGN §8.1), with molar stoichiometry ν
(`:H2 => -2, :O2 => -1, :H2O => 2`). `compile` checks that every reaction conserves every
element exactly. Extent rate ξ̇ [mol/s], affinity A = −Σ ν_i M_i μ_i [J/mol]. Bookkeeping:
m_i += ν_i M_i ξ̇, S += A ξ̇/T, production A ξ̇/T. So dU = 0 to round-off and the heat of
reaction appears through the formation enthalpies. `law` is an `AffinityLaw` (one for all
the reactions, coupled through L) or a `MassAction` (a vector of them for several reactions).
Connect it with `connect!(sys, reaction, body)`.
"""
struct Reaction{L}
    ν::Vector{Vector{Pair{Symbol,Float64}}}
    law::L
    function Reaction(νs::AbstractVector; law)
        isempty(νs) && throw(ArgumentError("a reaction element needs at least one reaction"))
        rs = Vector{Pair{Symbol,Float64}}[_stoichiometry(ν) for ν in νs]
        R = length(rs)
        if law isa AffinityLaw
            size(law.L, 1) == R || throw(ArgumentError("the affinity law is $(size(law.L, 1)) × $(size(law.L, 1)) for $R reactions"))
        elseif law isa MassAction
            R == 1 || throw(ArgumentError("give one MassAction per reaction (a vector) for $R reactions"))
            law = [law]
        elseif law isa AbstractVector && all(l -> l isa MassAction, law)
            length(law) == R || throw(ArgumentError("$(length(law)) mass-action laws for $R reactions"))
        else
            throw(ArgumentError("law must be an AffinityLaw, a MassAction or a vector of MassAction"))
        end
        return new{typeof(law)}(rs, law)
    end
end
Reaction(ν::Pair...; law) = Reaction([collect(ν)]; law)

function _stoichiometry(ν)
    out = Pair{Symbol,Float64}[Symbol(first(p)) => Float64(last(p)) for p in ν]
    isempty(out) && throw(ArgumentError("empty stoichiometry"))
    names = first.(out)
    allunique(names) || throw(ArgumentError("a species appears twice in one reaction: $names"))
    all(p -> isfinite(last(p)) && last(p) != 0, out) || throw(ArgumentError("stoichiometric coefficients must be finite and non-zero"))
    any(p -> last(p) < 0, out) && any(p -> last(p) > 0, out) ||
        throw(ArgumentError("a reaction needs reactants (ν < 0) and products (ν > 0)"))
    return out
end

# --- nuclear processes ------------------------------------------------------------------------------

"""
    FirstOrder(λ)
    FirstOrder(; half_life)

A first-order nuclear rate, ξ̇ = λ n_parent, where the parent is the single species with
ν = −1. λ [1/s] is a number or a function of time t -> λ(t) (σφ(t) for a neutron-driven rate).
"""
struct FirstOrder{F}
    λ::F
    function FirstOrder(λ)
        λ isa Real && !(λ >= 0 && isfinite(λ)) && throw(ArgumentError("decay constant must be ≥ 0 and finite, got $λ"))
        return new{typeof(λ)}(λ)
    end
end
function FirstOrder(; half_life::Real)
    (half_life > 0 && isfinite(half_life)) || throw(ArgumentError("half-life must be positive and finite"))
    return FirstOrder(log(2) / half_life)
end

"""
    NuclearSource(ν::Pair...; rate, E_dep = 0, external = (nucleons = 0, charge = 0))

A prescribed nuclear process in one body (decay, transmutation; DESIGN §8.1, D-11). `rate` is
`FirstOrder(λ)` or a function t -> ξ̇ [mol/s]. `E_dep` [J per mol of extent] is the energy
deposited in the body. `compile` checks Σ ν_i A_i = external.nucleons and Σ ν_i Z_i = external.charge, where
`external` declares what is absorbed from outside the model (the neutron of a capture:
nucleons = 1). Bookkeeping: m_i += ν_i M_i ξ̇ and dU = E_dep ξ̇ exactly, with
S += (E_dep ξ̇ + A ξ̇)/T. The species masses, the energy and the entropy are booked as external
(boundary) terms of either sign, never as production.
"""
struct NuclearSource{R}
    ν::Vector{Pair{Symbol,Float64}}
    rate::R
    E_dep::Float64
    external::NamedTuple{(:nucleons, :charge),Tuple{Int,Int}}
    parent::Symbol
    function NuclearSource(ν::Pair...; rate, E_dep::Real = 0.0, external = (nucleons = 0, charge = 0))
        st = _stoichiometry(ν)
        isfinite(E_dep) || throw(ArgumentError("E_dep must be finite"))
        ext = (nucleons = Int(get(external, :nucleons, 0)), charge = Int(get(external, :charge, 0)))
        parent = :none
        if rate isa FirstOrder
            ps = [first(p) for p in st if last(p) == -1]
            (length(ps) == 1 && count(p -> last(p) < 0, st) == 1) ||
                throw(ArgumentError("FirstOrder needs exactly one parent species, with ν = −1"))
            parent = ps[1]
        elseif !(rate isa Function)
            throw(ArgumentError("rate must be FirstOrder(λ) or a function t -> ξ̇"))
        end
        return new{typeof(rate)}(st, rate, Float64(E_dep), ext, parent)
    end
end

# --- heat sources -----------------------------------------------------------------------------------

"""
    HeatSource(Q̇)

A prescribed heat rate Q̇ [W] into a body (a number or a function t -> Q̇): S += Q̇/T. The energy
and the entropy Q̇/T are booked as external (boundary) terms, of either sign; it never touches
production. Nuclear (volumetric) heating, prescribed cooling.
"""
struct HeatSource{F}
    Q::F
    function HeatSource(Q)
        Q isa Real && !isfinite(Q) && throw(ArgumentError("heat rate must be finite, got $Q"))
        (Q isa Real || Q isa Function) || throw(ArgumentError("heat rate must be a number or a function of t"))
        return new{typeof(Q)}(Q)
    end
end
