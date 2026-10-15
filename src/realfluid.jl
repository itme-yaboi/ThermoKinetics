# Real fluids (DESIGN §2.4, §2.5; D-7, D-50). The backend lives in the package extension
# `ThermoKineticsClapeyronExt` (ext/), which loads with Clapeyron.jl; the core only declares the
# constructor and the hooks the rest of the library reads.

"""
    RealFluid(model; species)

A real-fluid material from a Clapeyron.jl equation of state (available once Clapeyron is
loaded; the backend is a package extension, D-50). Load it as `using ThermoKinetics;
import Clapeyron` and write `Clapeyron.IAPWS95()` etc.: Clapeyron also exports `pressure`,
`temperature` and `internal_energy`, which `using` both would make ambiguous; `import` costs
nothing at run time. `model` is any
Clapeyron `EoSModel` (`GERG2008(["helium"])`, `SingleFluid("carbon dioxide")`, `IAPWS95()`, a
mixture such as `GERG2008(["nitrogen", "oxygen", "argon"])`); `species` names its components
(a Symbol, or one Symbol per component, in the model's order).

Charges (S, V, m₁ … m_n), domains (Thermal, Volumetric, Species…), one coupling group, as for
an ideal gas. U(S, V, m) = A(T, V, n) + T S, with A = `Clapeyron.eos(model, V, T, n)` (the
volume-based Helmholtz energy, the only Clapeyron entry point of the state function) and T
from S = −∂A/∂T by safeguarded Newton; derivatives come from the implicit-function theorem, so
forces and the Hessian are exact (§2.4). A pure-fluid state inside the two-phase dome is a
saturated mixture (DESIGN §2.6): T from the lever rule, forces T, −p_sat, μ_sat, and `quality`
gives its vapour fraction; mixtures are single-phase only. The entropy and energy datums
are the model's own, so one species must not be held by two different backends in a model.

`charges_at(fluid; T, p or V, m)` builds charges (`m` a number for one component, else a
NamedTuple by species); `charges_at(fluid; T or p, x, m)` a saturated mixture of quality x.
"""
function RealFluid end

"""
The datum class of a material's species: materials exchanging a species must agree on it, or
the μ they exchange is not one force (D-48, D-50). `nothing` for a material without species.
"""
datum_class(::Material) = nothing
datum_class(::IdealGas) = :ideal_gas
datum_class(::SpeciesMaterial) = :species_data
datum_class(m::Moving) = datum_class(m.inner)

"Whether `isentropic_state(mat, s, p, y)` exists for a multi-species fluid (an ideal pump needs it)."
mixture_isentropic(::Material) = false

"`isentropic_state` with a composition: single-species materials ignore it."
isentropic_state(mat::Material, s, p, y) = isentropic_state(mat, s, p)

"""
    quality(material, E) -> Real or nothing

The vapour mass fraction of a two-phase state of a pure fluid (in [0, 1]), or `nothing` for a
single-phase state or a material without phases.
"""
quality(::Material, E) = nothing

"""
The specific entropy [J/(kg K)] of what a `draw` takes from a body with charges E: the mean
S/m for `:mixture` or a single-phase state, the saturated phase's otherwise (D-52).
"""
drawn_entropy(mat::Material, E, draw::Symbol) = E[thermal_slot(mat)] / sum(E[k] for k in mass_slots(mat))
