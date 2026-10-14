# Materials: one potential U(E) per material (DESIGN §2). Forces are its gradient and the
# equation-of-state matrix its Hessian, both by automatic differentiation, so P = ∂U/∂E and
# Maxwell symmetry hold by construction. Coupling groups are expressed by composition
# (`Composite`), never by flags.

abstract type Material end

"""
    domains(mat) -> Tuple

Ordered forms of motion of a body made of `mat`; its charges are stored in this order.
"""
function domains end

"""
    internal_energy(mat, E) -> U [J]

The material's state function, generic in the element type of `E` (for AD).
"""
function internal_energy end

ncharges(mat::Material) = sum(ncomponents, domains(mat))

"""
    forces(mat, E) -> SVector

P = ∂U/∂E. The only way any force is ever computed.
"""
forces(mat::Material, E::SVector) = ForwardDiff.gradient(e -> internal_energy(mat, e), E)
forces(mat::Material, E::AbstractVector) = forces(mat, SVector{ncharges(mat)}(E))

"""
    eos_matrix(mat, E) -> SMatrix

A = ∇²U, symmetrized (Maxwell reciprocity). For a material whose mass is a charge, A is
singular along E (first-order homogeneity); use `capacity_matrix` for capacities.
"""
function eos_matrix(mat::Material, E::AbstractVector)
    H = ForwardDiff.hessian(e -> internal_energy(mat, e), SVector{ncharges(mat)}(E))
    return (H + H') / 2
end

"Index of the thermal charge in `domains(mat)`, or `nothing` for a material without one."
thermal_slot(mat::Material) = findfirst(d -> d isa Thermal, domains(mat))

"Index of the volumetric charge in `domains(mat)`, or `nothing`."
volume_slot(mat::Material) = findfirst(d -> d isa Volumetric, domains(mat))

"Indices of the `Species` (mass) charges: the slots excluded from the fixed-mass block."
mass_slots(mat::Material) = Tuple(i for (i, d) in enumerate(domains(mat)) if d isa Species)

"Whether `E` is a state at which the material's forces are defined. Default: all finite."
valid_state(mat::Material, E::AbstractVector) = all(isfinite, E)

"""
    temperature(mat, E)

The thermal force T = ∂U/∂S.
"""
function temperature(mat::Material, E::AbstractVector)
    k = thermal_slot(mat)
    k === nothing && throw(ArgumentError("$(typeof(mat)) has no thermal form of motion"))
    return forces(mat, E)[k]
end

"""
    pressure(mat, E)

p = −∂U/∂V.
"""
function pressure(mat::Material, E::AbstractVector)
    k = volume_slot(mat)
    k === nothing && throw(ArgumentError("$(typeof(mat)) has no volumetric form of motion"))
    return -forces(mat, E)[k]
end

function _fixed_mass_block(mat::Material, E::AbstractVector)
    A = Matrix(eos_matrix(mat, E))
    idx = [i for i in 1:size(A, 1) if !(i in mass_slots(mat))]
    return A[idx, idx]
end

"""
    capacity_matrix(mat, E; hold = :forces) -> Matrix

Capacities K = dE/dP on the fixed-mass block (DESIGN §2.1, D-8; the full Hessian is never
inverted). `hold = :forces`: the other forces are held fixed (c_p-like), K = A_ff⁻¹.
`hold = :charges`: the other charges are held fixed (c_v-like), K_ii = 1/A_ii.
"""
function capacity_matrix(mat::Material, E::AbstractVector; hold::Symbol = :forces)
    A = _fixed_mass_block(mat, E)
    hold === :forces && return inv(A)
    hold === :charges && return Matrix(Diagonal(1 ./ diag(A)))
    throw(ArgumentError("hold must be :forces or :charges, got :$hold"))
end

"Thermodynamic stability: the fixed-mass block of ∇²U is positive definite."
is_stable(mat::Material, E::AbstractVector) = isposdef(Symmetric(_fixed_mass_block(mat, E)))

# --- constant-capacity solid ---------------------------------------------------------------

"""
    ConstantCapacitySolid(m, c; T_ref = 298.15)

An incompressible solid of fixed mass `m` [kg] and constant specific heat `c` [J/(kg K)].
One form of motion (thermal); the mass is a parameter, not a charge.

U(S) = m c T_ref exp(S / (m c)), so T = ∂U/∂S = T_ref exp(S / (m c)), U = m c T, and
S = 0 at T = T_ref.
"""
struct ConstantCapacitySolid <: Material
    m::Float64
    c::Float64
    T_ref::Float64
    function ConstantCapacitySolid(m::Real, c::Real; T_ref::Real = 298.15)
        (m > 0 && isfinite(m)) || throw(ArgumentError("mass must be positive and finite, got $m"))
        (c > 0 && isfinite(c)) || throw(ArgumentError("specific heat must be positive and finite, got $c"))
        (T_ref > 0 && isfinite(T_ref)) || throw(ArgumentError("T_ref must be positive and finite, got $T_ref"))
        return new(m, c, T_ref)
    end
end

domains(::ConstantCapacitySolid) = (Thermal(),)
internal_energy(s::ConstantCapacitySolid, E) = s.m * s.c * s.T_ref * exp(E[1] / (s.m * s.c))

"""
    charges_at(mat; kwargs...)

Charges of a body made of `mat` at a given state. Solids: `T`. Ideal gas: `T`, `m` and one
of `p` or `V`. Elastic volume: `V` or `p`.
"""
function charges_at(s::ConstantCapacitySolid; T::Real)
    _check_positive(T, "temperature")
    return SVector(s.m * s.c * log(T / s.T_ref))
end

function _check_positive(x, what)
    (x > 0 && isfinite(x)) || throw(ArgumentError("$what must be positive and finite, got $x"))
    return nothing
end

# --- ideal gas -----------------------------------------------------------------------------

"""
    IdealGas(species; R, cv, T_ref = 298.15, p_ref = 101325.0)

A single-species ideal gas with constant c_v: charges (S, V, m), domains
(Thermal, Volumetric, Species(species)). One coupled group (thermal + volumetric + mass).

Per unit mass, with v = V/m and s = S/m:
  T(s, v) = T_ref (v/v_ref)^(−R/c_v) exp(s/c_v),  v_ref = R T_ref / p_ref,
  U = m c_v T.
So s = 0 at (T_ref, p_ref), u = c_v T, h = c_p T and μ = h − T s. U is first-order
homogeneous in (S, V, m).
"""
struct IdealGas <: Material
    species::Symbol
    R::Float64
    cv::Float64
    T_ref::Float64
    p_ref::Float64
    function IdealGas(species::Symbol; R::Real, cv::Real, T_ref::Real = 298.15, p_ref::Real = 101325.0)
        _check_positive(R, "gas constant R")
        _check_positive(cv, "specific heat c_v")
        _check_positive(T_ref, "T_ref")
        _check_positive(p_ref, "p_ref")
        return new(species, R, cv, T_ref, p_ref)
    end
end

domains(g::IdealGas) = (Thermal(), Volumetric(), Species(g.species))

function internal_energy(g::IdealGas, E)
    S, V, m = E[1], E[2], E[3]
    v_ref = g.R * g.T_ref / g.p_ref
    T = g.T_ref * (V / (m * v_ref))^(-g.R / g.cv) * exp(S / (m * g.cv))
    return m * g.cv * T
end

valid_state(g::IdealGas, E::AbstractVector) = all(isfinite, E) && E[2] > 0 && E[3] > 0

function charges_at(g::IdealGas; T::Real, m::Real, p = nothing, V = nothing)
    _check_positive(T, "temperature")
    _check_positive(m, "mass")
    (p === nothing) == (V === nothing) && throw(ArgumentError("give exactly one of p or V"))
    if V === nothing
        _check_positive(p, "pressure")
        V = m * g.R * T / p
    else
        _check_positive(V, "volume")
    end
    v_ref = g.R * g.T_ref / g.p_ref
    s = g.cv * log(T / g.T_ref) + g.R * log((V / m) / v_ref)
    return SVector(m * s, Float64(V), Float64(m))
end

# --- elastic volume ------------------------------------------------------------------------

"""
    ElasticVolume(C; V0, p0 = 0.0)

An idealized elastic volume (a spring-loaded bladder or a stiff vessel wall) with one form of
motion, volumetric, and no thermal coupling: Veinik's linear one-DOF body P = A e (§11).
p = p0 + (V0 − V)/C, so the capacity dV/d(−p) is C > 0 [m³/Pa], and
U(V) = −p0 (V − V0) + (V − V0)²/(2C). It is **not** a gas: a gas must be modeled with its
thermal coupling (05, Boyle–Mariotte caution).
"""
struct ElasticVolume <: Material
    C::Float64
    V0::Float64
    p0::Float64
    function ElasticVolume(C::Real; V0::Real, p0::Real = 0.0)
        _check_positive(C, "capacity C")
        _check_positive(V0, "reference volume V0")
        isfinite(p0) || throw(ArgumentError("p0 must be finite"))
        return new(C, V0, p0)
    end
end

domains(::ElasticVolume) = (Volumetric(),)
internal_energy(e::ElasticVolume, E) = -e.p0 * (E[1] - e.V0) + (E[1] - e.V0)^2 / (2e.C)

function charges_at(e::ElasticVolume; V = nothing, p = nothing)
    (p === nothing) == (V === nothing) && throw(ArgumentError("give exactly one of p or V"))
    V === nothing && (V = e.V0 + e.C * (e.p0 - p))
    return SVector(Float64(V))
end

# --- composite (coupling groups) -------------------------------------------------------------

"""
    Composite(parts...)

A material made of independent parts on disjoint sets of domains: U = Σ U_part. Each part is
one coupling group, so ∇²U is block-diagonal by part, and a part's forces depend only on its
own charges (08 A7 holds by construction). At most one part may contain `Thermal`.
"""
struct Composite{P<:Tuple} <: Material
    parts::P
    function Composite(parts::Material...)
        isempty(parts) && throw(ArgumentError("a composite needs at least one part"))
        doms = FormOfMotion[d for p in parts for d in domains(p)]
        allunique(doms) || throw(ArgumentError("composite parts must have disjoint domains, got $doms"))
        count(d -> d isa Thermal, doms) <= 1 ||
            throw(ArgumentError("at most one part of a composite may contain the thermal form of motion"))
        return new{typeof(parts)}(parts)
    end
end

domains(c::Composite) = _cat_domains(c.parts)
_cat_domains(::Tuple{}) = ()
_cat_domains(parts::Tuple) = (domains(first(parts))..., _cat_domains(Base.tail(parts))...)

@inline _part_slice(E, off, n) = SVector{n}(ntuple(k -> E[off+k], n))

internal_energy(c::Composite, E) = _energy_parts(c.parts, E, 0)
_energy_parts(::Tuple{}, E, off) = zero(eltype(E))
function _energy_parts(parts::Tuple, E, off)
    p = first(parts)
    n = ncharges(p)
    return internal_energy(p, _part_slice(E, off, n)) + _energy_parts(Base.tail(parts), E, off + n)
end

valid_state(c::Composite, E::AbstractVector) = _valid_parts(c.parts, E, 0)
_valid_parts(::Tuple{}, E, off) = true
function _valid_parts(parts::Tuple, E, off)
    p = first(parts)
    n = ncharges(p)
    return valid_state(p, _part_slice(E, off, n)) && _valid_parts(Base.tail(parts), E, off + n)
end

"""
    charges_at(c::Composite, states::NamedTuple...)

One NamedTuple of `charges_at` keywords per part, in order, e.g.
`charges_at(Composite(solid, spring), (T = 300.0,), (p = 2e5,))`.
"""
function charges_at(c::Composite, states::NamedTuple...)
    length(states) == length(c.parts) ||
        throw(ArgumentError("give one state per part ($(length(c.parts))), got $(length(states))"))
    return vcat((charges_at(p; s...) for (p, s) in zip(c.parts, states))...)
end
