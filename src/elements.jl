# Bodies, reservoirs and conductors (DESIGN §3, §5). Elements hold parameters and laws only;
# they never write state derivatives (that is the assembler's job, assembly.jl).

"""
    Body(name, material, E₀)
    Body(name, material; T)

A body: a material plus its initial charges. Its charges are its state.
"""
struct Body{M<:Material,N}
    name::Symbol
    material::M
    E0::SVector{N,Float64}
end

function Body(name::Symbol, mat::Material, E0::AbstractVector)
    N = ncharges(mat)
    length(E0) == N || throw(ArgumentError("$(typeof(mat)) has $N charges, got $(length(E0))"))
    return Body{typeof(mat),N}(name, mat, SVector{N,Float64}(E0))
end
Body(name::Symbol, mat::Material; T::Real) = Body(name, mat, charges_at(mat; T))

"""
    Reservoir(name; T)

Fixed forces, infinite capacity. In slice 1 only the thermal force is supported.
"""
struct Reservoir
    name::Symbol
    T::Float64
    function Reservoir(name::Symbol; T::Real)
        (T > 0 && isfinite(T)) || throw(ArgumentError("reservoir temperature must be positive and finite, got $T"))
        return new(name, T)
    end
end

has_thermal(b::Body) = thermal_slot(b.material) !== nothing
has_thermal(::Reservoir) = true

# --- conductor laws ---------------------------------------------------------------------

abstract type ConductorLaw end

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

"Kernel A coefficient M_QQ for the thermal row at end temperatures `Ta`, `Tb`."
conductance_QQ(law::ThermalConductance, Ta, Tb) = law.G * Ta * Tb

domains(::ThermalConductance) = (Thermal(),)

"""
    Conductor(law)

A conductor element; joined to two nodes with `connect!`.
"""
struct Conductor{L<:ConductorLaw}
    law::L
end

domains(c::Conductor) = domains(c.law)
