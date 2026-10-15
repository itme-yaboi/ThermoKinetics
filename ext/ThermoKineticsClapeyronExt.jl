"""
    ThermoKineticsClapeyronExt

The real-fluid backend (DESIGN §2.4, §2.5; D-7, D-50): `RealFluid(model; species)` for any
Clapeyron.jl equation of state. Loaded automatically with `using ThermoKinetics, Clapeyron`.

The state function uses only `Clapeyron.eos(model, V, T, n)`, the volume-based Helmholtz
energy, and its automatic derivatives (noted by the user, D-50). T is found from the entropy
by safeguarded Newton on Float64 values; dual-number inputs (forces, Hessians) get their
partials from one Newton step per nesting level taken in their own number type, which is the
implicit-function theorem: exact derivatives, with AD never running through the iterations.
"""
module ThermoKineticsClapeyronExt

using ThermoKinetics
using Clapeyron
using ForwardDiff
using StaticArrays
import ThermoKinetics: RealFluid, domains, internal_energy, valid_state, charges_at, isentropic_state,
                       datum_class, mixture_isentropic, species_names, molar_mass, Material

const R̄ = 8.314462618

"A Clapeyron-backed fluid: the model, its components' names and molar masses [kg/mol]."
struct ClapeyronFluid{Mod,N} <: Material
    model::Mod
    names::NTuple{N,Symbol}
    M::NTuple{N,Float64}
    Tc::Float64                        # critical temperature (pure fluids), NaN for mixtures
    Tlast::Base.RefValue{Float64}      # warm start for the inversion
end

function RealFluid(model::Clapeyron.EoSModel; species)
    names = species isa Symbol ? (species,) : Tuple(Symbol.(species))
    M = Tuple(Float64.(Clapeyron.mw(model)) ./ 1000)
    length(names) == length(M) ||
        throw(ArgumentError("the model has $(length(M)) components, but $(length(names)) species names were given"))
    allunique(names) || throw(ArgumentError("species names must be distinct"))
    Tc = length(M) == 1 ? Float64(first(Clapeyron.crit_pure(model))) : NaN
    return ClapeyronFluid{typeof(model),length(M)}(model, names, M, Tc, Ref(NaN))
end

Base.show(io::IO, f::ClapeyronFluid) = print(io, "RealFluid(", nameof(typeof(f.model)), "; species = ", f.names, ")")

domains(f::ClapeyronFluid) = (Thermal(), Volumetric(), map(Species, f.names)...)
species_names(f::ClapeyronFluid) = f.names
molar_mass(f::ClapeyronFluid, name::Symbol) = f.M[findfirst(==(name), f.names)]
datum_class(f::ClapeyronFluid) = (:clapeyron, nameof(typeof(f.model)))
mixture_isentropic(::ClapeyronFluid) = true

# --- the Helmholtz energy and its temperature derivatives (the only Clapeyron entry point) ------

@inline _moles(f::ClapeyronFluid{Mod,N}, E) where {Mod,N} = ntuple(k -> E[2+k] / f.M[k], Val(N))
_helm(f::ClapeyronFluid, V, T, n) = Clapeyron.eos(f.model, V, T, collect(promote(n...)))
_entropy(f, T, V, n) = -ForwardDiff.derivative(t -> _helm(f, V, t, n), T)
_dSdT(f, T, V, n) = -ForwardDiff.derivative(t -> ForwardDiff.derivative(t2 -> _helm(f, V, t2, n), t), T)
_pressure(f, T, V, n) = -ForwardDiff.derivative(v -> _helm(f, v, T, n), V)

"Safeguarded Newton for S(T; V, n) = S on Float64 values: bracket by halving or doubling, then Newton with a bisection fallback."
function _solve_T(f::ClapeyronFluid, S::Float64, V::Float64, n::NTuple{N,Float64}) where {N}
    T = f.Tlast[]
    (isfinite(T) && T > 0) || (T = 300.0)
    g(T) = _entropy(f, T, V, n) - S
    ok(T) = (d = _dSdT(f, T, V, n); isfinite(d) && d > 0 && isfinite(g(T)))   # inside the model's monotone range
    ok(T) || (T = 300.0)
    lo, hi = T, T
    # bracket by steps of 1.2, staying where ∂S/∂T > 0: below its range (e.g. under the triple point)
    # an equation of state extrapolates wildly and S(T) is no longer monotone
    # (a step that lands outside that range is shrunk, not taken: near the critical point a full step
    # from above it can jump into the dome)
    if g(T) > 0
        found, fac = false, 1.2
        for _ in 1:400
            Tn = lo / fac
            if !ok(Tn)
                fac = sqrt(fac)
                fac < 1 + 1e-9 ? break : continue
            end
            lo = Tn
            g(lo) < 0 && (found = true; break)
        end
        found || throw(DomainError(S, "no temperature with this entropy at V = $V within the model's range (below it)"))
    else
        found, fac = false, 1.2
        for _ in 1:400
            Tn = hi * fac
            if !ok(Tn)
                fac = sqrt(fac)
                fac < 1 + 1e-9 ? break : continue
            end
            hi = Tn
            g(hi) > 0 && (found = true; break)
        end
        found || throw(DomainError(S, "no temperature with this entropy at V = $V within the model's range (above it)"))
    end
    T = (lo + hi) / 2
    for _ in 1:200
        gT = g(T)
        gT > 0 ? (hi = T) : (lo = T)
        d = _dSdT(f, T, V, n)
        Tn = T - gT / d
        (isfinite(Tn) && lo < Tn < hi && d > 0) || (Tn = (lo + hi) / 2)
        if abs(Tn - T) <= 4eps() * T
            T = Tn
            break
        end
        T = Tn
    end
    f.Tlast[] = T
    return T
end

"T(S, V, n) for any number types: solve on the values, then one Newton step in the inputs' type (IFT)."
function _temperature(f::ClapeyronFluid, S, V, n::Tuple)
    Tp = promote_type(typeof(S), typeof(V), map(typeof, n)...)
    Tp === Float64 && return _solve_T(f, Float64(S), Float64(V), map(Float64, n))
    T0 = _temperature(f, ForwardDiff.value(S), ForwardDiff.value(V), map(ForwardDiff.value, n))
    return T0 - (_entropy(f, T0, V, n) - S) / _dSdT(f, T0, V, n)
end

function internal_energy(f::ClapeyronFluid{Mod,N}, E) where {Mod,N}
    S, V = E[1], E[2]
    n = _moles(f, E)
    T = _temperature(f, S, V, n)
    return _helm(f, V, T, n) + T * S
end

"The Float64 value of a number, through any nesting of dual numbers (validity needs values only)."
_val(x::ForwardDiff.Dual) = _val(ForwardDiff.value(x))
_val(x::Real) = Float64(x)

function valid_state(f::ClapeyronFluid{Mod,N}, Ein::AbstractVector) where {Mod,N}
    E = map(_val, Ein)
    all(isfinite, E) && E[2] > 0 || return false
    all(k -> E[2+k] >= 0, 1:N) && sum(E[3:end]) > 0 || return false
    n = _moles(f, E)
    T = try
        _temperature(f, E[1], E[2], n)
    catch err
        err isa DomainError || rethrow()
        return false
    end
    # a pure fluid inside the two-phase dome is not a single-phase state (slice 10b brings the flash)
    if N == 1 && T < f.Tc
        sat = Clapeyron.saturation_pressure(f.model, T)
        v = E[2] / sum(n)
        (isfinite(sat[2]) && isfinite(sat[3]) && sat[2] < v < sat[3]) && return false
    end
    return true
end

function _masses(f::ClapeyronFluid{Mod,N}, m) where {Mod,N}
    ms = if m isa Real
        N == 1 || throw(ArgumentError("give the masses of a $N-component fluid as a NamedTuple by species"))
        (Float64(m),)
    else
        Set(keys(m)) == Set(f.names) || throw(ArgumentError("amounts must name exactly the species $(f.names)"))
        ntuple(k -> Float64(m[f.names[k]]), N)
    end
    all(x -> x >= 0 && isfinite(x), ms) && sum(ms) > 0 || throw(ArgumentError("masses must be ≥ 0, finite, with a positive total"))
    return ms
end

"""
    charges_at(fluid::RealFluid; T, p or V, m)

Charges at temperature T and either pressure p or volume V. With p, the volume is Clapeyron's
stable root `volume(model, p, T, n)` (a constructor convenience: the state function itself uses
only the volume-based Helmholtz energy).
"""
function charges_at(f::ClapeyronFluid{Mod,N}; T::Real, p = nothing, V = nothing, m) where {Mod,N}
    T > 0 || throw(ArgumentError("temperature must be positive"))
    count(!isnothing, (p, V)) == 1 || throw(ArgumentError("give T, m and exactly one of p, V"))
    ms = _masses(f, m)
    n = ntuple(k -> ms[k] / f.M[k], N)
    if V === nothing
        p > 0 || throw(ArgumentError("pressure must be positive"))
        V = _volume_pT(f, Float64(p), Float64(T), n)
    end
    V > 0 || throw(ArgumentError("volume must be positive"))
    S = _entropy(f, Float64(T), Float64(V), n)
    return SVector{N + 2,Float64}(S, V, ms...)
end

"""
Clapeyron's stable-root volume at (p, T) for amounts n, computed for one mole of the composition
and scaled: its root finder returns NaN for compressed liquids at small amounts (~1e-4 mol).
"""
function _volume_pT(f::ClapeyronFluid, p::Float64, T::Float64, n::Tuple)
    ntot = sum(n)
    return ntot * Clapeyron.volume(f.model, p, T, collect(n) ./ ntot)
end

# --- the isentropic state at a pressure (pumps, turbines; §2.4 item 8) ---------------------------

"""
    isentropic_state(fluid, s, p[, y]) -> (; T, μ)

The state of one kilogram (composition `y`, mass fractions) at specific entropy s and pressure
p: a 2-D Newton in (ln T, ln V) on S(T, V) = s and p(T, V) = p using only the Helmholtz energy,
from a guess at Clapeyron's (p, T) root; dual inputs by one Newton step in their type (IFT).
μ is per kg, a number for one component and a tuple otherwise.
"""
isentropic_state(f::ClapeyronFluid{Mod,1}, s, p) where {Mod} = isentropic_state(f, s, p, (1.0,))

function isentropic_state(f::ClapeyronFluid{Mod,N}, s, p, y) where {Mod,N}
    n = ntuple(k -> y[k] / f.M[k], Val(N))                   # moles in one kilogram
    T, V = _isentropic_TV(f, s, p, n)
    E = SVector{N + 2}(s, V, ntuple(k -> y[k] * one(V), Val(N))...)
    F = ThermoKinetics.forces(f, E)
    μ = ntuple(k -> F[2+k], Val(N))
    return (; T = F[1], μ = N == 1 ? μ[1] : μ)
end

function _isentropic_TV(f::ClapeyronFluid, s, p, n::Tuple)
    Tp = promote_type(typeof(s), typeof(p), map(typeof, n)...)
    if Tp === Float64
        return _solve_isentropic(f, Float64(s), Float64(p), map(Float64, n))
    end
    T0, V0 = _isentropic_TV(f, ForwardDiff.value(s), ForwardDiff.value(p), map(ForwardDiff.value, n))
    # one Newton step on r(T, V) = (S − s, p − p_target) in the inputs' number type
    r1 = _entropy(f, T0, V0, n) - s
    r2 = _pressure(f, T0, V0, n) - p
    J = ForwardDiff.jacobian(x -> [_entropy(f, x[1], x[2], n), _pressure(f, x[1], x[2], n)], [T0, V0])
    d = J \ [r1, r2]
    return T0 - d[1], V0 - d[2]
end

function _solve_isentropic(f::ClapeyronFluid, s::Float64, p::Float64, n::NTuple{N,Float64}) where {N}
    # 1. locate the branch: a bracketed 1-D solve in T along the isobar, with Clapeyron's stable root
    #    V(p, T) (S increases with T along an isobar); 2. polish on the volume-based equations alone
    Vp(T) = _volume_pT(f, p, T, n)
    h(T) = _entropy(f, T, Vp(T), n) - s
    T = isfinite(f.Tlast[]) ? f.Tlast[] : 300.0
    lo, hi = T, T
    if h(T) > 0
        for _ in 1:200
            lo /= 1.2
            (isfinite(h(lo)) && h(lo) < 0) && break
            lo < 1.0 && throw(DomainError(s, "no isentropic state at p = $p within the model's range"))
        end
    else
        for _ in 1:200
            hi *= 1.2
            (isfinite(h(hi)) && h(hi) > 0) && break
            hi > 1e5 && throw(DomainError(s, "no isentropic state at p = $p within the model's range"))
        end
    end
    for _ in 1:100                                          # bisection to ~1e-10 relative
        mid = (lo + hi) / 2
        h(mid) > 0 ? (hi = mid) : (lo = mid)
        (hi - lo) <= 1e-10 * hi && break
    end
    T = (lo + hi) / 2
    x = [log(T), log(Vp(T))]
    for _ in 1:100
        T, V = exp(x[1]), exp(x[2])
        r = [_entropy(f, T, V, n) - s, (_pressure(f, T, V, n) - p) / p]
        J = ForwardDiff.jacobian(z -> [_entropy(f, exp(z[1]), exp(z[2]), n), _pressure(f, exp(z[1]), exp(z[2]), n) / p], x)
        dx = -(J \ r)
        λ = 1.0                                              # damping: at most a factor e per step in T or V
        while maximum(abs, λ * dx) > 1.0
            λ /= 2
        end
        x .+= λ .* dx
        maximum(abs, λ * dx) <= 8eps() && break
    end
    T, V = exp(x[1]), exp(x[2])
    abs(_entropy(f, T, V, n) - s) <= 1e-9 * abs(s) + 1e-9 && abs(_pressure(f, T, V, n) - p) <= 1e-9 * p ||
        throw(DomainError(s, "no single-phase isentropic state at p = $p (it may lie in the two-phase region)"))
    return T, V
end

end # module
