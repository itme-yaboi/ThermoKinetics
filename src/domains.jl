# Forms of motion: declarations with no physics (DESIGN §1). Physics lives in material
# potentials (static) and conductor laws (kinetic). There are no sign flags: every force
# is P = ∂U/∂E.

abstract type FormOfMotion end

"Thermal form of motion: charge entropy S [J/K], force temperature T [K]."
struct Thermal <: FormOfMotion end

"""
Volumetric (fluid-mechanical) form of motion: charge volume V [m³], force −p [Pa].
The force is −p because P = ∂U/∂V = −p; there is no sign exception anywhere.
"""
struct Volumetric <: FormOfMotion end

"""
    Species(name)

Chemical form of motion for one species: charge mass m [kg], force chemical potential
μ [J/kg]. Two `Species` with the same name are the same domain.
"""
struct Species <: FormOfMotion
    name::Symbol
end

charge_symbol(::Thermal) = :S
charge_unit(::Thermal) = "J/K"
force_symbol(::Thermal) = :T
force_unit(::Thermal) = "K"

charge_symbol(::Volumetric) = :V
charge_unit(::Volumetric) = "m^3"
force_symbol(::Volumetric) = Symbol("-p")
force_unit(::Volumetric) = "Pa"

charge_symbol(::Species) = :m
charge_unit(::Species) = "kg"
force_symbol(::Species) = :μ
force_unit(::Species) = "J/kg"

"Number of state slots the charge occupies (vector charges such as momentum take more)."
ncomponents(::FormOfMotion) = 1

"Whether the charge can be produced. Only thermal charge (entropy) can."
producible(::FormOfMotion) = false
producible(::Thermal) = true

"""
Default absolute force tolerance used to scale integrator tolerances: each charge is
resolved to the equivalent of this force change, through its capacity.
"""
default_force_atol(::Thermal) = 1e-8      # K
default_force_atol(::Volumetric) = 1e-4   # Pa
default_force_atol(::Species) = 1e-6      # J/kg
