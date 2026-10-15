# Moving bodies of variable mass (DESIGN §19.2, D-33) and the frozen-gas properties that pumps
# and the nozzle read from a fluid material.

"""
    Moving(material; m_dry = 0.0)

A body that moves along one axis while holding mass: the inner `material`'s charges followed by
a `Momentum` charge p, with

  U(E, p) = U_int(E) + p²/(2M),  M = m_dry + Σ_k m_k,

so the force on p is the velocity w = p/M, T and p are those of the inner material (Galilean
invariance), and μ_k = μ_k,int − w²/2. Mass and momentum are one coupling group. Every body that
exchanges mass with a moving body must move with it (a stream or pump between a moving and a
non-moving body is refused at compile); bodies of one vehicle are held at one velocity by
`Transformer(1; a = Momentum(), b = Momentum())`. U is first-order homogeneous only for
m_dry = 0. Build a body with `Body(:tank, Moving(gas; m_dry); w = 0.0, state...)`.
"""
struct Moving{M<:Material,NI,NS} <: Material
    inner::M
    m_dry::Float64
    ms::NTuple{NS,Int}          # the inner material's mass slots
    function Moving(inner::Material; m_dry::Real = 0.0)
        any(d -> d isa Momentum, domains(inner)) && throw(ArgumentError("the material already has a Momentum charge"))
        ms = mass_slots(inner)
        isempty(ms) && throw(ArgumentError("Moving needs a material with mass (Species charges)"))
        (m_dry >= 0 && isfinite(m_dry)) || throw(ArgumentError("dry mass must be ≥ 0 and finite, got $m_dry"))
        return new{typeof(inner),ncharges(inner),length(ms)}(inner, Float64(m_dry), ms)
    end
end

domains(m::Moving) = (domains(m.inner)..., Momentum())

@inline _inner_charges(m::Moving{M,NI}, E) where {M,NI} = SVector{NI}(ntuple(k -> E[k], Val(NI)))
@inline _total_mass(m::Moving{M,NI,NS}, E) where {M,NI,NS} = m.m_dry + sum(ntuple(k -> E[m.ms[k]], Val(NS)))

function internal_energy(m::Moving{M,NI}, E) where {M,NI}
    p = E[NI+1]
    return internal_energy(m.inner, _inner_charges(m, E)) + p^2 / (2 * _total_mass(m, E))
end

valid_state(m::Moving{M,NI}, E::AbstractVector) where {M,NI} =
    isfinite(E[NI+1]) && valid_state(m.inner, view(E, 1:NI)) && _total_mass(m, E) > 0

"""
    charges_at(mv::Moving; w = 0.0, state...)

The inner material's charges at `state...`, followed by the momentum M w.
"""
function charges_at(m::Moving; w::Real = 0.0, state...)
    isfinite(w) || throw(ArgumentError("velocity must be finite, got $w"))
    Ei = charges_at(m.inner; state...)
    M = m.m_dry + sum(Ei[k] for k in m.ms)
    return vcat(Ei, SVector(M * w))
end

stiffness_floor(m::Moving) = (stiffness_floor(m.inner)..., 0.0)
species_names(m::Moving) = species_names(m.inner)
species_data(m::Moving, name::Symbol) = species_data(m.inner, name)
molar_mass(m::Moving, name::Symbol) = molar_mass(m.inner, name)

"The trace floor of a material's mixing term (D-41), or `nothing` for a material without one."
trace_floor(::Material) = nothing
trace_floor(m::SpeciesMaterial) = m.trace
trace_floor(m::Moving) = trace_floor(m.inner)

"The material that holds the fluid's thermodynamics (the inner one of a `Moving` body)."
fluid_material(mat::Material) = mat
fluid_material(m::Moving) = m.inner

# --- isentropic states and frozen-gas properties ------------------------------------------------

isentropic_state(m::Moving, s, p) = isentropic_state(m.inner, s, p)

"""
    isentropic_state(mix::IdealGasMixture, s, p) -> (; T, μ)

For a one-species mixture: the temperature and chemical potential at specific entropy `s` and
pressure `p`. For one species the material's entropy is exactly s = s0 + c_p ln(T/T_ref) −
R ln(p/p_ref) (the species is the whole mass, far above the trace floor), so T is closed-form;
μ then comes from the material itself, at one kilogram.
"""
function isentropic_state(g::IdealGasMixture{1}, s, p)
    cp = g.cv[1] + g.R[1]
    T = g.T_ref * exp((s - g.s0[1] + g.R[1] * log(p / g.p_ref)) / cp)
    F = forces(g, SVector(s, g.R[1] * T / p, one(s)))
    return (; T = F[1], μ = F[3])
end

"""
    frozen_gas(mat, E) -> (; R, cv, cp, k, m)

Mass-weighted gas constant and specific heats of a gas body at frozen composition [J/(kg K)], and
its fluid mass. Constant per species, so the frozen expansion is calorically perfect.
"""
frozen_gas(g::IdealGas, E) = (; R = g.R, cv = g.cv, cp = g.cv + g.R, k = (g.cv + g.R) / g.cv, m = E[3])
function frozen_gas(g::IdealGasMixture{N}, E) where {N}
    m = sum(ntuple(k -> E[2+k], Val(N)))
    R = sum(ntuple(k -> E[2+k] * g.R[k], Val(N))) / m
    cv = sum(ntuple(k -> E[2+k] * g.cv[k], Val(N))) / m
    return (; R, cv, cp = cv + R, k = (cv + R) / cv, m)
end
frozen_gas(m::Moving, E) = frozen_gas(m.inner, _inner_charges(m, E))
frozen_gas(mat::Material, E) = throw(ArgumentError("a nozzle needs an ideal-gas chamber; $(typeof(mat)) has no frozen-gas properties"))
