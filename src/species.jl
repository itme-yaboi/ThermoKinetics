# Species data and multi-species materials (DESIGN §8.1): ideal-gas mixtures and ideal
# condensed solutions with formation enthalpies, so a reaction's affinity carries its heat.

"Molar gas constant R̄ [J/(mol K)]."
const R_MOLAR = 8.314462618

"""
Atomic data per symbol: molar mass [kg/mol] and nucleon number. Element symbols use the IUPAC
standard atomic weights and the nucleon number of their main isotope; isotope symbols (D, T,
He3, …) use their atomic masses. `n` is the free neutron. Species outside this table give
their own `atom_data`.
"""
const ATOMIC_DATA = Dict{Symbol,NamedTuple{(:M, :A),Tuple{Float64,Int}}}(
    :H => (M = 1.008e-3, A = 1), :He => (M = 4.0026e-3, A = 4), :Li => (M = 6.94e-3, A = 7),
    :Be => (M = 9.0122e-3, A = 9), :B => (M = 10.81e-3, A = 11), :C => (M = 12.011e-3, A = 12),
    :N => (M = 14.007e-3, A = 14), :O => (M = 15.999e-3, A = 16), :F => (M = 18.998e-3, A = 19),
    :Ne => (M = 20.180e-3, A = 20), :Na => (M = 22.990e-3, A = 23), :Mg => (M = 24.305e-3, A = 24),
    :Al => (M = 26.982e-3, A = 27), :Si => (M = 28.085e-3, A = 28), :P => (M = 30.974e-3, A = 31),
    :S => (M = 32.06e-3, A = 32), :Cl => (M = 35.45e-3, A = 35), :Ar => (M = 39.95e-3, A = 40),
    :K => (M = 39.098e-3, A = 39), :Ca => (M = 40.078e-3, A = 40), :Ti => (M = 47.867e-3, A = 48),
    :Cr => (M = 51.996e-3, A = 52), :Mn => (M = 54.938e-3, A = 55), :Fe => (M = 55.845e-3, A = 56),
    :Co => (M = 58.933e-3, A = 59), :Ni => (M = 58.693e-3, A = 58), :Cu => (M = 63.546e-3, A = 63),
    :Zr => (M = 91.224e-3, A = 90), :Mo => (M = 95.95e-3, A = 98), :Ag => (M = 107.87e-3, A = 107),
    :I => (M = 126.90e-3, A = 127), :Xe => (M = 131.29e-3, A = 132), :W => (M = 183.84e-3, A = 184),
    :Pb => (M = 207.2e-3, A = 208), :U => (M = 238.03e-3, A = 238),
    # isotopes and the neutron
    :n => (M = 1.00866491595e-3, A = 1), :D => (M = 2.01410177812e-3, A = 2), :T => (M = 3.01604928e-3, A = 3),
    :He3 => (M = 3.01602932e-3, A = 3), :He4 => (M = 4.00260325413e-3, A = 4),
    :Li6 => (M = 6.0151228874e-3, A = 6), :Li7 => (M = 7.0160034366e-3, A = 7),
    :B10 => (M = 10.01293695e-3, A = 10), :B11 => (M = 11.00930536e-3, A = 11),
)

"""
    SpeciesData(name; atoms, cv, h_f = 0, s0 = 0, charge = 0, basis = :mass, atom_data = (;))

One species: its atoms (`(H = 2, O = 1)`; isotopes are their own symbols), its constant specific
heat `cv` (at constant volume for a gas; the specific heat of a condensed phase), formation
enthalpy `h_f` at T_ref and standard entropy `s0` at (T_ref, p_ref). With `basis = :mass` these
are per kg [J/(kg K)], [J/kg]; with `basis = :molar`, per mole, and they are stored per kg.
The molar mass and nucleon number come from the atoms (`ATOMIC_DATA`, or `atom_data =
(X = (M = …, A = …),)` for symbols outside it), so reactions conserve mass by construction.
`charge` is the net electric charge (0 for a neutral species). Data are phase-specific: liquid
and gaseous water are two `SpeciesData` with the same name and atoms.
"""
struct SpeciesData{A<:NamedTuple}
    name::Symbol
    atoms::A                  # symbol => count
    atom_M::Dict{Symbol,Float64}
    M::Float64                # kg/mol
    nucleons::Int
    charge::Int
    cv::Float64               # J/(kg K)
    h_f::Float64              # J/kg
    s0::Float64               # J/(kg K)
end

function SpeciesData(name::Symbol; atoms::NamedTuple, cv::Real, h_f::Real = 0.0, s0::Real = 0.0, charge::Integer = 0,
                     basis::Symbol = :mass, atom_data::NamedTuple = (;))
    isempty(atoms) && throw(ArgumentError("species :$name needs at least one atom"))
    M = 0.0
    A = 0
    atom_M = Dict{Symbol,Float64}()
    for (sym, cnt) in pairs(atoms)
        (cnt isa Integer && cnt > 0) || throw(ArgumentError("species :$name: atom counts must be positive integers, got $sym = $cnt"))
        d = haskey(atom_data, sym) ? atom_data[sym] : get(ATOMIC_DATA, sym, nothing)
        d === nothing && throw(ArgumentError("species :$name: unknown atom :$sym; give its data with atom_data = ($sym = (M = …, A = …),)"))
        (d.M > 0 && isfinite(d.M)) || throw(ArgumentError("species :$name: atom :$sym needs a positive molar mass"))
        atom_M[sym] = Float64(d.M)
        M += cnt * d.M
        A += cnt * d.A
    end
    (cv > 0 && isfinite(cv)) || throw(ArgumentError("species :$name: specific heat must be positive and finite, got $cv"))
    (isfinite(h_f) && isfinite(s0)) || throw(ArgumentError("species :$name: h_f and s0 must be finite"))
    basis in (:mass, :molar) || throw(ArgumentError("basis must be :mass or :molar, got :$basis"))
    f = basis === :molar ? 1 / M : 1.0
    return SpeciesData{typeof(atoms)}(name, atoms, atom_M, M, A, Int(charge), cv * f, h_f * f, s0 * f)
end

"Molar mass [kg/mol] of a species."
molar_mass(sd::SpeciesData) = sd.M
"Nucleon number of a species."
nucleons(sd::SpeciesData) = sd.nucleons
"Atom count of `sym` in a species (0 if absent)."
atom_count(sd::SpeciesData, sym::Symbol) = haskey(sd.atoms, sym) ? Int(sd.atoms[sym]) : 0

# --- shared helpers for multi-species materials -----------------------------------------------

"Materials whose species carry `SpeciesData` (reactions and nuclear sources need it)."
abstract type SpeciesMaterial <: Material end

species_names(mat::SpeciesMaterial) = mat.names
function species_data(mat::SpeciesMaterial, name::Symbol)
    k = findfirst(==(name), mat.names)
    k === nothing && throw(ArgumentError("no species :$name in this material (it holds $(mat.names))"))
    return mat.species[k]
end
species_data(::Material, name::Symbol) = nothing
molar_mass(mat::SpeciesMaterial, name::Symbol) = molar_mass(species_data(mat, name))
species_names(::Material) = ()

function _species_tuple(species)
    isempty(species) && throw(ArgumentError("a mixture needs at least one species"))
    names = Tuple(sd.name for sd in species)
    allunique(names) || throw(ArgumentError("duplicate species in a mixture: $names"))
    return names
end

"Masses from `m` or moles `n`, as NamedTuples keyed by species name (any order)."
function _species_masses(mat::SpeciesMaterial, m, n)
    count(!isnothing, (m, n)) == 1 || throw(ArgumentError("give the species amounts as either m (kg) or n (mol)"))
    amt = m === nothing ? n : m
    Set(keys(amt)) == Set(mat.names) ||
        throw(ArgumentError("amounts must name exactly the species $(mat.names), got $(keys(amt))"))
    out = ntuple(k -> Float64(amt[mat.names[k]]) * (m === nothing ? mat.M[k] : 1.0), length(mat.names))
    all(x -> x >= 0 && isfinite(x), out) || throw(ArgumentError("species amounts must be ≥ 0 and finite"))
    sum(out) > 0 || throw(ArgumentError("a mixture needs a positive total amount"))
    return out
end


# --- the trace floor (D-41) ----------------------------------------------------------------------
#
# Ideal mixing has n ln n, whose derivative ln n → −∞ as n → 0. A species whose equilibrium
# amount lies below round-off (a burnt-out fuel, a nuclide after many half-lives) would then sit
# in round-off noise at the edge of the valid states, and steps that land slightly below zero
# would be rejected without end. Below a floor δ (`trace` × the body's total mass), n ln n is
# continued by its second-order Taylor expansion at δ: U stays a convex, C², first-order
# homogeneous potential, μ is finite for every amount, and a depleted species settles at a
# small finite value (it can read slightly negative) instead of at a singular edge. Species
# above the floor are exactly ideal.

"n ln n, continued below the floor δ by its Taylor expansion at δ."
@inline function _xlogx(n, δ)
    n >= δ && return n * log(n)
    d = n - δ
    return δ * log(δ) + (log(δ) + 1) * d + d^2 / (2δ)
end

"ln of the activity amount, consistent with `_xlogx`: ln n above the floor, ln δ + (n − δ)/δ below."
@inline _log_activity(n, δ) = n >= δ ? log(n) : log(δ) + (n - δ) / δ

@inline _total_mass(E, off, ::Val{N}) where {N} = sum(ntuple(k -> E[off+k], Val(N)))

function _check_trace(trace)
    (trace > 0 && trace < 1e-6) || throw(ArgumentError("trace floor must be in (0, 1e-6), got $trace"))
    return Float64(trace)
end

"A state is valid if its total mass is positive and no species is below −1e-6 of it (a sanity bound)."
@inline function _valid_amounts(E, off, ::Val{N}) where {N}
    mt = _total_mass(E, off, Val(N))
    return mt > 0 && all(k -> E[off+k] > -1e-6 * mt, 1:N)
end

# --- ideal-gas mixture --------------------------------------------------------------------------

"""
    IdealGasMixture(species::SpeciesData...; T_ref = 298.15, p_ref = 101325.0, trace = 1e-13)

An ideal-gas mixture with constant c_v per species: charges (S, V, m₁ … m_n), domains
(Thermal, Volumetric, Species(name)…), one coupling group. With R_i = R̄/M_i,
v_ref,i = R_i T_ref/p_ref and C_v = Σ m_i c_v,i:
  ln(T/T_ref) = (S − Σ m_i [s0_i + R_i ln(V/(m_i v_ref,i))]) / C_v,
  U = Σ m_i [h_f,i − R_i T_ref + c_v,i (T − T_ref)].
So p = Σ m_i R_i T/V (Dalton) and μ_i = h_i(T) − T s_i(T, p_i), the mixing entropy included.
U is first-order homogeneous.

A species may be absent or depleted. Below `trace` × Σ m (default 1e-13 of the body's mass),
the mixing term n ln n is continued smoothly (D-41), so μ stays finite and a depleted species
can read slightly negative instead of failing. Above it the mixture is exactly ideal.
It matters only if a species below the floor has a driving force you need (its permeation,
a rate that depends on it, its equilibrium pressure): then lower `trace`, but not below
~1e-15 (round-off), or give that species its own body (D-42).
"""
struct IdealGasMixture{N,S<:Tuple} <: SpeciesMaterial
    species::S
    names::NTuple{N,Symbol}
    M::NTuple{N,Float64}
    R::NTuple{N,Float64}
    cv::NTuple{N,Float64}
    hf::NTuple{N,Float64}
    s0::NTuple{N,Float64}
    vref::NTuple{N,Float64}
    T_ref::Float64
    p_ref::Float64
    trace::Float64
    function IdealGasMixture(species::SpeciesData...; T_ref::Real = 298.15, p_ref::Real = 101325.0, trace::Real = 1e-13)
        names = _species_tuple(species)
        _check_positive(T_ref, "T_ref")
        _check_positive(p_ref, "p_ref")
        N = length(species)
        M = ntuple(k -> species[k].M, N)
        R = ntuple(k -> R_MOLAR / M[k], N)
        return new{N,typeof(species)}(species, names, M, R, ntuple(k -> species[k].cv, N), ntuple(k -> species[k].h_f, N),
                                      ntuple(k -> species[k].s0, N), ntuple(k -> R[k] * T_ref / p_ref, N), T_ref, p_ref,
                                      _check_trace(trace))
    end
end

domains(g::IdealGasMixture) = (Thermal(), Volumetric(), map(Species, g.names)...)

"Σ m_i s0_i + R̄ Σ [n_i ln(V/(M_i v_ref,i)) − n_i ln n_i] (the configurational entropy at T_ref)."
function _gas_mixture_σ(g::IdealGasMixture{N}, V, m) where {N}
    δ = g.trace * sum(m)
    return sum(ntuple(k -> m[k] * g.s0[k] +
                           R_MOLAR * ((m[k] / g.M[k]) * log(V / (g.M[k] * g.vref[k])) - _xlogx(m[k] / g.M[k], δ / g.M[k])),
                      Val(N)))
end

function internal_energy(g::IdealGasMixture{N}, E) where {N}
    m = ntuple(k -> E[2+k], Val(N))
    Cv = sum(ntuple(k -> m[k] * g.cv[k], Val(N)))
    T = g.T_ref * exp((E[1] - _gas_mixture_σ(g, E[2], m)) / Cv)
    return sum(ntuple(k -> m[k] * (g.hf[k] - g.R[k] * g.T_ref + g.cv[k] * (T - g.T_ref)), Val(N)))
end

valid_state(g::IdealGasMixture{N}, E::AbstractVector) where {N} = all(isfinite, E) && E[2] > 0 && _valid_amounts(E, 2, Val(N))

"""
    charges_at(mix::IdealGasMixture; T, V or p, m or n)

Charges at temperature T with species amounts `m` (kg) or `n` (mol), as NamedTuples by
species name (zero allowed), and either the volume V or the total pressure p.
"""
function charges_at(g::IdealGasMixture{N}; T::Real, V = nothing, p = nothing, m = nothing, n = nothing) where {N}
    _check_positive(T, "temperature")
    ms = _species_masses(g, m, n)
    count(!isnothing, (V, p)) == 1 || throw(ArgumentError("give T, the amounts and exactly one of V, p"))
    if V === nothing
        _check_positive(p, "pressure")
        V = sum(ms[k] * g.R[k] for k in 1:N) * T / p
    end
    _check_positive(V, "volume")
    S = sum(ms[k] * g.cv[k] for k in 1:N) * log(T / g.T_ref) + _gas_mixture_σ(g, Float64(V), ms)
    return SVector{N + 2,Float64}(S, V, ms...)
end

"""
    log_activity_amounts(mat, E) -> NTuple

ln of the amounts [mol] in the mixing terms: ln n_i above the trace floor, its continuation
below (D-41). Mass action uses them, so its rates stay consistent with μ.
"""
function log_activity_amounts(g::IdealGasMixture{N}, E) where {N}
    δ = g.trace * _total_mass(E, 2, Val(N))
    return ntuple(k -> _log_activity(E[2+k] / g.M[k], δ / g.M[k]), Val(N))
end

# --- ideal condensed solution -------------------------------------------------------------------

"""
    IdealSolution(species::SpeciesData...; T_ref = 298.15, trace = 1e-13)

A rigid condensed mixture (a liquid or a solid in a vessel) with ideal mixing on mole
fractions: charges (S, m₁ … m_n), domains (Thermal, Species(name)…), no `Volumetric` (it is
incompressible and its pressure is not a state). With c_i the species' specific heat, C = Σ m_i c_i
and x_i the mole fractions:
  ln(T/T_ref) = (S − Σ m_i [s0_i − R_i ln x_i]) / C,
  U = Σ m_i [h_f,i + c_i (T − T_ref)],
so μ_i = h_i(T) − T s_i(T, x_i), molar μ°_i + R̄T ln x_i. U is first-order homogeneous in (S, m).
Below the trace floor the mixing term is continued smoothly, as for `IdealGasMixture` (D-41).
"""
struct IdealSolution{N,S<:Tuple} <: SpeciesMaterial
    species::S
    names::NTuple{N,Symbol}
    M::NTuple{N,Float64}
    R::NTuple{N,Float64}
    cv::NTuple{N,Float64}
    hf::NTuple{N,Float64}
    s0::NTuple{N,Float64}
    T_ref::Float64
    trace::Float64
    function IdealSolution(species::SpeciesData...; T_ref::Real = 298.15, trace::Real = 1e-13)
        names = _species_tuple(species)
        _check_positive(T_ref, "T_ref")
        N = length(species)
        M = ntuple(k -> species[k].M, N)
        return new{N,typeof(species)}(species, names, M, ntuple(k -> R_MOLAR / M[k], N), ntuple(k -> species[k].cv, N),
                                      ntuple(k -> species[k].h_f, N), ntuple(k -> species[k].s0, N), Float64(T_ref),
                                      _check_trace(trace))
    end
end

domains(s::IdealSolution) = (Thermal(), map(Species, s.names)...)

"Σ m_i s0_i − R̄ [Σ n_i ln n_i − n ln n]: the standard and mixing entropy at T_ref."
function _solution_σ(s::IdealSolution{N}, m) where {N}
    δ = s.trace * sum(m)
    ntot = sum(ntuple(k -> m[k] / s.M[k], Val(N)))
    mix = sum(ntuple(k -> _xlogx(m[k] / s.M[k], δ / s.M[k]), Val(N))) - ntot * log(ntot)
    return sum(ntuple(k -> m[k] * s.s0[k], Val(N))) - R_MOLAR * mix
end

function internal_energy(s::IdealSolution{N}, E) where {N}
    m = ntuple(k -> E[1+k], Val(N))
    C = sum(ntuple(k -> m[k] * s.cv[k], Val(N)))
    T = s.T_ref * exp((E[1] - _solution_σ(s, m)) / C)
    return sum(ntuple(k -> m[k] * (s.hf[k] + s.cv[k] * (T - s.T_ref)), Val(N)))
end

valid_state(s::IdealSolution{N}, E::AbstractVector) where {N} = all(isfinite, E) && _valid_amounts(E, 1, Val(N))

"""
    charges_at(sol::IdealSolution; T, m or n)

Charges at temperature T with species amounts `m` (kg) or `n` (mol), as NamedTuples by name
(zero allowed).
"""
function charges_at(s::IdealSolution{N}; T::Real, m = nothing, n = nothing) where {N}
    _check_positive(T, "temperature")
    ms = _species_masses(s, m, n)
    S = sum(ms[k] * s.cv[k] for k in 1:N) * log(T / s.T_ref) + _solution_σ(s, ms)
    return SVector{N + 1,Float64}(S, ms...)
end

function log_activity_amounts(s::IdealSolution{N}, E) where {N}
    δ = s.trace * _total_mass(E, 1, Val(N))
    return ntuple(k -> _log_activity(E[1+k] / s.M[k], δ / s.M[k]), Val(N))
end
