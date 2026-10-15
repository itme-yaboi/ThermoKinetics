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

"""
Electrical form of motion: charge q [C], force electric potential φ [V]. Across material
boundaries the true force for carriers is the electrochemical potential (`04`); within one
conductor material φ serves. φ has no absolute zero, so only differences enter any law.
"""
struct Electrical <: FormOfMotion end

"""
Kinetic translation (one axis): charge momentum p [kg m/s], force velocity w [m/s]
(Veinik (7); 04). The inertial view (inertia force × displacement, (9)) is the same
phenomenon seen from outside, so it is never a separate domain.
"""
struct Momentum <: FormOfMotion end

"Kinetic rotation: charge angular momentum L [kg m²/s], force angular velocity ω [rad/s] ((8))."
struct AngularMomentum <: FormOfMotion end

"""
The kinetic form of the electrical domain: charge flux linkage Λ [Wb = V s], force current
i [A]. An inductor stores it, and an inductor branch converts it to charge (a 3-port converter).
"""
struct FluxLinkage <: FormOfMotion end

charge_symbol(::Thermal) = :S
charge_unit(::Thermal) = "J/K"
force_symbol(::Thermal) = :T
force_unit(::Thermal) = "K"

charge_symbol(::Volumetric) = :V
charge_unit(::Volumetric) = "m^3"
force_symbol(::Volumetric) = Symbol("-p")
force_unit(::Volumetric) = "Pa"

charge_symbol(::Electrical) = :q
charge_unit(::Electrical) = "C"
force_symbol(::Electrical) = :φ
force_unit(::Electrical) = "V"

charge_symbol(::Species) = :m
charge_unit(::Species) = "kg"
force_symbol(::Species) = :μ
force_unit(::Species) = "J/kg"

charge_symbol(::Momentum) = :p
charge_unit(::Momentum) = "kg m/s"
force_symbol(::Momentum) = :w
force_unit(::Momentum) = "m/s"

charge_symbol(::AngularMomentum) = :L
charge_unit(::AngularMomentum) = "kg m^2/s"
force_symbol(::AngularMomentum) = :ω
force_unit(::AngularMomentum) = "rad/s"

charge_symbol(::FluxLinkage) = :Λ
charge_unit(::FluxLinkage) = "Wb"
force_symbol(::FluxLinkage) = :i
force_unit(::FluxLinkage) = "A"

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
default_force_atol(::Electrical) = 1e-8   # V
default_force_atol(::Momentum) = 1e-9     # m/s
default_force_atol(::AngularMomentum) = 1e-8   # rad/s
default_force_atol(::FluxLinkage) = 1e-9  # A
