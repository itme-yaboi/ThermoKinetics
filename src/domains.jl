# Forms of motion: declarations with no physics (DESIGN §1). Physics lives in material
# potentials (static) and conductor laws (kinetic). There are no sign flags: every force
# is P = ∂U/∂E.

abstract type FormOfMotion end

"Thermal form of motion: charge entropy S [J/K], force temperature T [K]."
struct Thermal <: FormOfMotion end

charge_symbol(::Thermal) = :S
charge_unit(::Thermal) = "J/K"
force_symbol(::Thermal) = :T
force_unit(::Thermal) = "K"

"Number of state slots the charge occupies (vector charges such as momentum take more)."
ncomponents(::FormOfMotion) = 1

"Whether the charge can be produced. Only thermal charge (entropy) can."
producible(::FormOfMotion) = false
producible(::Thermal) = true
