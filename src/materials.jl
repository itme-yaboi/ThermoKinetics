# Materials: one potential U(E) per material (DESIGN §2). Forces are its gradient and the
# equation-of-state matrix its Hessian, both by automatic differentiation, so P = ∂U/∂E and
# Maxwell symmetry hold by construction.

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

A = ∇²U, symmetric by construction (Maxwell reciprocity).
"""
function eos_matrix(mat::Material, E::AbstractVector)
    H = ForwardDiff.hessian(e -> internal_energy(mat, e), SVector{ncharges(mat)}(E))
    return (H + H') / 2
end

"Index of the thermal charge in `domains(mat)`, or `nothing` for a body without one."
thermal_slot(mat::Material) = findfirst(d -> d isa Thermal, domains(mat))

"""
    temperature(mat, E)

The thermal force T = ∂U/∂S.
"""
function temperature(mat::Material, E::AbstractVector)
    k = thermal_slot(mat)
    k === nothing && throw(ArgumentError("$(typeof(mat)) has no thermal form of motion"))
    return forces(mat, E)[k]
end

# ---------------------------------------------------------------------------------------

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
    charges_at(mat; T)

Charges of a body made of `mat` at temperature `T` [K].
"""
function charges_at(s::ConstantCapacitySolid; T::Real)
    (T > 0 && isfinite(T)) || throw(ArgumentError("temperature must be positive and finite, got $T"))
    return SVector(s.m * s.c * log(T / s.T_ref))
end
