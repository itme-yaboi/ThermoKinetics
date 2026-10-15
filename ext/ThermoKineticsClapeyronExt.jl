"""
    ThermoKineticsClapeyronExt

The real-fluid backend (DESIGN §2.4, §2.5; D-7, D-50): `RealFluid(model; species)` for any
Clapeyron.jl equation of state. Loaded automatically with `using ThermoKinetics, Clapeyron`.

The state function uses only `Clapeyron.eos(model, V, T, n)`, the volume-based Helmholtz
energy, and its automatic derivatives (noted by the user, D-50). T is found from the entropy
by safeguarded Newton on Float64 values; dual-number inputs (forces, Hessians) get their
partials from one Newton step per nesting level taken in their own number type, which is the
implicit-function theorem: exact derivatives, with AD never running through the iterations.

A pure fluid inside its two-phase dome is a saturated mixture (slice 10b, DESIGN §2.6): its
temperature comes from the lever rule, its saturated volumes from equal pressure and equal μ
on the same Helmholtz energy, and derivatives again from one Newton step per dual level.
"""
module ThermoKineticsClapeyronExt

using ThermoKinetics
using Clapeyron
using ForwardDiff
using StaticArrays
using LinearAlgebra: Symmetric, eigvals
import ThermoKinetics: RealFluid, domains, internal_energy, valid_state, charges_at, isentropic_state,
                       datum_class, mixture_isentropic, species_names, molar_mass, Material,
                       quality, drawn_entropy, is_stable

const R̄ = 8.314462618

"A Clapeyron-backed fluid: the model, its components' names and molar masses [kg/mol]."
struct ClapeyronFluid{Mod,N} <: Material
    model::Mod
    names::NTuple{N,Symbol}
    M::NTuple{N,Float64}
    Tc::Float64                        # critical temperature (pure fluids), NaN for mixtures
    pc::Float64                        # critical pressure (pure fluids), NaN for mixtures
    Tlo::Float64                       # low end of the saturation range (pure fluids), NaN for mixtures
    plo::Float64                       # the saturation pressure there
    Tlast::Base.RefValue{Float64}      # warm start for the inversion
    T2last::Base.RefValue{Float64}     # warm start for the flash
    satlast::Base.RefValue{NTuple{3,Float64}}    # (T, v_l, v_v) of the last saturation solve
    Tsatlast::Base.RefValue{NTuple{2,Float64}}   # (p, T_sat) of the last
end

function RealFluid(model::Clapeyron.EoSModel; species)
    names = species isa Symbol ? (species,) : Tuple(Symbol.(species))
    M = Tuple(Float64.(Clapeyron.mw(model)) ./ 1000)
    length(names) == length(M) ||
        throw(ArgumentError("the model has $(length(M)) components, but $(length(names)) species names were given"))
    allunique(names) || throw(ArgumentError("species names must be distinct"))
    Tc, pc, Tlo, plo = NaN, NaN, NaN, NaN
    if length(M) == 1
        cp = Clapeyron.crit_pure(model)
        Tc, pc = Float64(cp[1]), Float64(cp[2])
        # the saturation range starts at 0.4 T_c (below the triple point of most fluids) or where the
        # model's saturation solver first succeeds
        for fr in 0.4:0.02:0.9
            ps = first(Clapeyron.saturation_pressure(model, fr * Tc))
            if isfinite(ps)
                Tlo, plo = fr * Tc, Float64(ps)
                break
            end
        end
    end
    return ClapeyronFluid{typeof(model),length(M)}(model, names, M, Tc, pc, Tlo, plo, Ref(NaN), Ref(NaN),
                                                    Ref((NaN, NaN, NaN)), Ref((NaN, NaN)))
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

"""
T(S, V, n) for any number types: solve on the values (or take `Tval`, already solved there),
then one Newton step in the inputs' type per dual level (IFT).
"""
function _temperature(f::ClapeyronFluid, S, V, n::Tuple, Tval::Float64 = NaN)
    Tp = promote_type(typeof(S), typeof(V), map(typeof, n)...)
    Tp === Float64 && return isfinite(Tval) ? Tval : _solve_T(f, Float64(S), Float64(V), map(Float64, n))
    T0 = _temperature(f, ForwardDiff.value(S), ForwardDiff.value(V), map(ForwardDiff.value, n), Tval)
    return T0 - (_entropy(f, T0, V, n) - S) / _dSdT(f, T0, V, n)
end

function internal_energy(f::ClapeyronFluid{Mod,N}, E) where {Mod,N}
    S, V = E[1], E[2]
    n = _moles(f, E)
    Tval = NaN
    if N == 1                                    # a pure fluid may be a saturated mixture
        two, Tval = _phase(f, _val(S), _val(V), (_val(n[1]),))
        two && return _u2(f, S, V, n, Tval)
    end
    T = _temperature(f, S, V, n, Tval)
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
    try
        N == 1 ? _phase(f, E[1], E[2], n) : _temperature(f, E[1], E[2], n)
    catch err
        err isa DomainError || rethrow()
        return false
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
    charges_at(fluid::RealFluid; T or p, x, m)

Charges at temperature T and either pressure p or volume V. With p, the volume is Clapeyron's
stable root `volume(model, p, T, n)` (a constructor convenience: the state function itself uses
only the volume-based Helmholtz energy). With a quality `x` (a pure fluid), the saturated
mixture at T or at the saturation temperature of p.
"""
function charges_at(f::ClapeyronFluid{Mod,N}; T = nothing, p = nothing, V = nothing, x = nothing, m) where {Mod,N}
    x === nothing || return _charges_wet(f, T, p, V, x, m)
    T === nothing && throw(ArgumentError("give T, m and exactly one of p, V (or T or p with a quality x)"))
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

function _charges_wet(f::ClapeyronFluid{Mod,N}, T, p, V, x, m) where {Mod,N}
    N == 1 || throw(ArgumentError("a quality x is defined for a pure fluid only"))
    (V === nothing && count(!isnothing, (T, p)) == 1) ||
        throw(ArgumentError("with a quality x, give exactly one of T, p (and m)"))
    0 <= x <= 1 || throw(ArgumentError("the quality must lie in [0, 1]"))
    Ts = T === nothing ? _Tsat_values(f, Float64(p)) : Float64(T)
    r = Ts === nothing ? nothing : _sat_values(f, Ts)
    r === nothing && throw(ArgumentError("no saturation state at T = $T, p = $p within the model's range"))
    ms = _masses(f, m)
    n = ms[1] / f.M[1]
    sl, sv = _sat_entropies(f, Ts, r[1], r[2])
    return SVector{3,Float64}(n * ((1 - x) * sl + x * sv), n * ((1 - x) * r[1] + x * r[2]), ms[1])
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
μ is per kg, a number for one component and a tuple otherwise. For a pure fluid with
s_l(p) < s < s_v(p) the state is wet: T_sat(p) and the saturated μ (a turbine's wet exhaust).
"""
isentropic_state(f::ClapeyronFluid{Mod,1}, s, p) where {Mod} = isentropic_state(f, s, p, (1.0,))

function isentropic_state(f::ClapeyronFluid{Mod,N}, s, p, y) where {Mod,N}
    if N == 1
        w = _wet_isentropic(f, s, p)
        w === nothing || return w
    end
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

# --- two-phase states of a pure fluid (slice 10b; DESIGN §2.6, D-52) ----------------------------
# Everything below works per mole (n = 1) on the same volume-based Helmholtz energy.

const _ONE = (1.0,)

"The saturation residual at T in z = (ln v_l, ln v_v): equal pressures, equal chemical potentials."
function _sat_residual(f::ClapeyronFluid, z, T)
    vl, vv = exp(z[1]), exp(z[2])
    al, av = _helm(f, vl, T, _ONE), _helm(f, vv, T, _ONE)
    pl, pv = _pressure(f, T, vl, _ONE), _pressure(f, T, vv, _ONE)
    RT = R̄ * T
    return SVector((pl - pv) * vv / RT, (al + pl * vl - av - pv * vv) / RT)
end

_sat_newton(f::ClapeyronFluid, z, T) = z - ForwardDiff.jacobian(w -> _sat_residual(f, w, T), z) \ _sat_residual(f, z, T)

"Newton on the saturation residual from z; the converged z, or nothing (no convergence, or one root for both phases)."
function _sat_polish(f::ClapeyronFluid, z::SVector{2,Float64}, T::Float64)
    for _ in 1:12
        zn = _sat_newton(f, z, T)
        all(isfinite, zn) || return nothing
        d = maximum(abs, zn - z)
        z = zn
        if d <= 1e-13
            # two distinct, mechanically stable phases (∂p/∂v < 0): not the trivial root, nor a
            # spurious one in the unstable region
            z[2] - z[1] > 1e-6 || return nothing
            dpdv(v) = ForwardDiff.derivative(w -> _pressure(f, T, w, _ONE), v)
            (dpdv(exp(z[1])) < 0 && dpdv(exp(z[2])) < 0) || return nothing
            return z
        end
    end
    return nothing
end

"""
Saturated molar volumes (v_l, v_v) at T (Float64), or `nothing` outside [T_lo, T_c). Newton from
the last solve when T is within 0.1 % of it, otherwise from Clapeyron's `saturation_pressure`;
either way the result must stay near its start (a multiparameter equation of state has spurious
equal-p, equal-μ pairs at high density, which a distant start can fall into).
"""
function _sat_values(f::ClapeyronFluid, T::Float64)
    c = f.satlast[]
    c[1] == T && return (c[2], c[3])
    (f.Tlo <= T < f.Tc) || return nothing
    z = nothing
    near(z, z0, tol) = z !== nothing && maximum(abs, z - z0) <= tol
    if isfinite(c[1]) && abs(T - c[1]) <= 1e-3 * T && max(T, c[1]) < 0.99 * f.Tc
        z0 = SVector(log(c[2]), log(c[3]))
        z = _sat_polish(f, z0, T)
        near(z, z0, 0.05) || (z = nothing)
    end
    if z === nothing
        _, vl, vv = Clapeyron.saturation_pressure(f.model, T)
        (isfinite(vl) && isfinite(vv) && 0 < vl < vv) || return nothing
        z0 = SVector(log(Float64(vl)), log(Float64(vv)))
        z = _sat_polish(f, z0, T)
        near(z, z0, 1e-3) || return nothing
    end
    vl, vv = exp(z[1]), exp(z[2])
    f.satlast[] = (T, vl, vv)
    return (vl, vv)
end

"(v_l, v_v) at T for any number type: the values, then one Newton step per dual level (IFT)."
function _sat(f::ClapeyronFluid, T)
    T isa ForwardDiff.Dual || return _sat_values(f, Float64(T))
    r = _sat(f, ForwardDiff.value(T))
    r === nothing && return nothing
    z = _sat_newton(f, SVector(log(r[1]), log(r[2])), T)
    return (exp(z[1]), exp(z[2]))
end

"Molar entropies (s_l, s_v) of the saturated phases at T."
_sat_entropies(f::ClapeyronFluid, T, vl, vv) = (_entropy(f, T, vl, _ONE), _entropy(f, T, vv, _ONE))

"The saturation pressure at T, or nothing."
_psat(f::ClapeyronFluid, T) = (r = _sat(f, T); r === nothing ? nothing : _pressure(f, T, r[1], _ONE))

"""
The lever-rule residual at T for molar (s, v): the two-phase entropy on the isochore v minus s
(it increases with T). `nothing` where v is outside (v_l, v_v) (above the isochore's dome
temperature) unless `check = false`.
"""
function _lever(f::ClapeyronFluid, T, s, v; check::Bool = true)
    r = _sat(f, T)
    r === nothing && return nothing
    vl, vv = r
    check && !(vl < v < vv) && return nothing
    sl, sv = _sat_entropies(f, T, vl, vv)
    return sl + (v - vl) / (vv - vl) * (sv - sl) - s
end

"""
The two-phase temperature of molar (s, v) (Float64): safeguarded Newton on the lever rule in
[T_lo, T_c]. A DomainError when (s, v) is not inside the dome within that range.
"""
function _flash_values(f::ClapeyronFluid, s::Float64, v::Float64)
    lo, hi = f.Tlo, f.Tc
    g = _lever(f, lo, s, v)
    (g !== nothing && g < 0) || throw(DomainError((s, v), "not a two-phase state within the saturation range"))
    T = (isfinite(f.T2last[]) && lo < f.T2last[] < hi) ? f.T2last[] : (lo + hi) / 2
    inside = false
    for _ in 1:200
        g = _lever(f, T, s, v)
        if g === nothing                         # v outside the dome at T: the root lies below
            hi = T
            T = (lo + hi) / 2
            continue
        end
        inside = true
        g > 0 ? (hi = T) : (lo = T)
        d = ForwardDiff.derivative(t -> _lever(f, t, s, v; check = false), T)
        Tn = T - g / d
        (isfinite(Tn) && lo < Tn < hi && d > 0) || (Tn = (lo + hi) / 2)
        if abs(Tn - T) <= 4eps() * T || hi - lo <= 4eps() * T
            T = Tn
            break
        end
        T = Tn
    end
    # a root, not the isochore's dome boundary (where the bracket also closes, with g < 0, for a
    # state that is single-phase)
    r = inside ? _sat_values(f, T) : nothing
    ok = r !== nothing && -1e-9 <= (v - r[1]) / (r[2] - r[1]) <= 1 + 1e-9 &&
         abs(_lever(f, T, s, v; check = false)) <= 1e-8 * max(abs(s), R̄)
    ok || throw(DomainError((s, v), "not a two-phase state within the saturation range"))
    f.T2last[] = T
    return T
end

"""
The phase of a pure-fluid state (Float64): `(false, T)` single phase, `(true, T)` two phase.
The single-phase inversion comes first. Its state is the stable one unless its volume lies
between the saturated volumes at its own temperature (a metastable or unstable state); then,
or when no single-phase temperature exists, the flash decides. A DomainError if neither applies.
"""
function _phase(f::ClapeyronFluid{Mod,1}, S::Float64, V::Float64, n::Tuple{Float64}) where {Mod}
    v = V / n[1]
    T1 = try
        _solve_T(f, S, V, n)
    catch err
        err isa DomainError || rethrow()
        NaN
    end
    if isfinite(T1)
        T1 >= f.Tc && return (false, T1)
        r = _sat_values(f, T1)
        if r === nothing
            # T1 below the saturation range (a metastable root of the inversion) or at the critical
            # end: the state is single-phase unless the flash finds it in the dome
            try
                return (true, _flash_values(f, S / n[1], v))
            catch err
                err isa DomainError || rethrow()
                return (false, T1)
            end
        end
        (v <= r[1] || v >= r[2]) && return (false, T1)
    end
    return (true, _flash_values(f, S / n[1], v))
end

"The two-phase T for any number types: `Tval` on the values, then one Newton step per dual level (IFT)."
function _T2(f::ClapeyronFluid, S, V, n::Tuple, Tval::Float64)
    promote_type(typeof(S), typeof(V), typeof(n[1])) <: ForwardDiff.Dual || return Tval
    T0 = _T2(f, ForwardDiff.value(S), ForwardDiff.value(V), map(ForwardDiff.value, n), Tval)
    s, v = S / n[1], V / n[1]
    return T0 - _lever(f, T0, s, v; check = false) / ForwardDiff.derivative(t -> _lever(f, t, s, v; check = false), T0)
end

"U of a saturated mixture: A + T S with A = n[(1 − x) a(v_l, T) + x a(v_v, T)], x by the lever rule."
function _u2(f::ClapeyronFluid, S, V, n::Tuple, Tval::Float64)
    T = _T2(f, S, V, n, Tval)
    vl, vv = _sat(f, T)
    x = (V / n[1] - vl) / (vv - vl)
    return n[1] * ((1 - x) * _helm(f, vl, T, _ONE) + x * _helm(f, vv, T, _ONE)) + T * S
end

function quality(f::ClapeyronFluid{Mod,1}, Ein::AbstractVector) where {Mod}
    E = map(_val, Ein)
    n = _moles(f, E)
    two, T = _phase(f, E[1], E[2], n)
    two || return nothing
    vl, vv = _sat_values(f, T)
    return clamp((E[2] / n[1] - vl) / (vv - vl), 0.0, 1.0)
end

function drawn_entropy(f::ClapeyronFluid{Mod,1}, E, draw::Symbol) where {Mod}
    S, V, m = E[1], E[2], E[3]
    draw === :mixture && return S / m
    n = _moles(f, E)
    two, Tval = _phase(f, _val(S), _val(V), (_val(n[1]),))
    two || return S / m
    T = _T2(f, S, V, n, Tval)
    vl, vv = _sat(f, T)
    return _entropy(f, T, draw === :liquid ? vl : vv, _ONE) / f.M[1]
end

"""
Inside the dome the fixed-mass Hessian is singular along the tie line (the phase rule): the
state is stable when that block is positive semidefinite.
"""
function is_stable(f::ClapeyronFluid{Mod,1}, Ein::AbstractVector) where {Mod}
    E = map(_val, Ein)
    two = try
        first(_phase(f, E[1], E[2], _moles(f, E)))
    catch err
        err isa DomainError || rethrow()
        false
    end
    two || return invoke(is_stable, Tuple{Material,AbstractVector}, f, Ein)
    λ = eigvals(Symmetric(Matrix(ThermoKinetics._fixed_mass_block(f, E))))
    return minimum(λ) >= -1e-8 * maximum(abs, λ)
end

"""
T_sat(p) (Float64) by safeguarded Newton on ln p_sat(T) = ln p in [T_lo, T_c]; `nothing` when
p is outside the saturation range.
"""
function _Tsat_values(f::ClapeyronFluid, p::Float64)
    c = f.Tsatlast[]
    c[1] == p && return c[2]
    lo, hi = f.Tlo, f.Tc
    (isfinite(lo) && f.plo < p < f.pc) || return nothing
    # the start: the last solve, or ln p linear in 1/T between (T_lo, p_lo) and the critical point
    T = (isfinite(c[2]) && lo < c[2] < hi) ? c[2] :
        1 / (1 / lo + (log(p) - log(f.plo)) / (log(f.pc) - log(f.plo)) * (1 / hi - 1 / lo))
    for _ in 1:200
        ps = _psat(f, T)
        if ps === nothing                        # at the critical end, where the phases merge
            hi = T
            T = (lo + hi) / 2
            continue
        end
        g = log(ps) - log(p)
        g > 0 ? (hi = T) : (lo = T)
        d = ForwardDiff.derivative(t -> log(_psat(f, t)), T)
        Tn = T - g / d
        (isfinite(Tn) && lo < Tn < hi && d > 0) || (Tn = (lo + hi) / 2)
        if abs(Tn - T) <= 4eps() * T || hi - lo <= 4eps() * T
            T = Tn
            break
        end
        T = Tn
    end
    f.Tsatlast[] = (p, T)
    return T
end

"T_sat(p) for any number type: `Tval` on the values, then one Newton step per dual level (IFT)."
function _Tsat(f::ClapeyronFluid, p, Tval::Float64)
    p isa ForwardDiff.Dual || return Tval
    T0 = _Tsat(f, ForwardDiff.value(p), Tval)
    return T0 - (log(_psat(f, T0)) - log(p)) / ForwardDiff.derivative(t -> log(_psat(f, t)), T0)
end

"""
The isentropic state of one kilogram of a pure fluid at (s, p) when it is wet,
s_l(p) < s < s_v(p): T_sat(p) and the saturated μ per kg. `nothing` otherwise.
"""
function _wet_isentropic(f::ClapeyronFluid{Mod,1}, s, p) where {Mod}
    Tv = _Tsat_values(f, _val(p))
    Tv === nothing && return nothing
    r = _sat_values(f, Tv)
    r === nothing && return nothing
    sl, sv = _sat_entropies(f, Tv, r[1], r[2]) ./ f.M[1]
    sl < _val(s) < sv || return nothing
    T = _Tsat(f, p, Tv)
    vl, _ = _sat(f, T)
    μ = (_helm(f, vl, T, _ONE) + _pressure(f, T, vl, _ONE) * vl) / f.M[1]
    z = zero(s) * zero(p)                        # the common number type of the inputs
    return (; T = T + z, μ = μ + z)
end

end # module
