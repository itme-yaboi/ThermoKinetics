# Spatial helpers (`07` §6): a 1D finite-volume chain is bodies (cells) plus
# conductors (faces). Nothing here is special to the solver.

"""
    rod!(sys, name; n, L, area, λ, ρ, c, T0, left = nothing, right = nothing, T_ref = 298.15)

Add a uniform 1D rod of `n` equal cells (constant-capacity solid) to `sys`. `T0` is a number
or a function of the cell-centre position x [m]. Each boundary is one of:
- `nothing`: insulated;
- a reservoir label: fixed temperature at the rod end (half-cell conductance 2λA/Δx);
- `(reservoir_label, h)`: convection with coefficient h [W/(m² K)] in series with the
  half cell, G = 1/(1/(hA) + Δx/(2λA)).

Positive heat flow runs from left to right. Returns
`(; cells, faces, left_face, right_face, all_faces, x, dx)`, where `all_faces` lists every
conductor in order from left to right.
"""
function rod!(sys::ThermoSystem, name::Symbol; n::Integer, L::Real, area::Real, λ::Real, ρ::Real, c::Real,
              T0, left = nothing, right = nothing, T_ref::Real = 298.15)
    n >= 1 || throw(ArgumentError("a rod needs at least one cell"))
    dx = L / n
    x = [(i - 0.5) * dx for i in 1:n]
    mat = ConstantCapacitySolid(ρ * area * dx, c; T_ref)
    Tinit(xi) = T0 isa Function ? T0(xi) : T0
    cells = [add!(sys, Body(Symbol(name, :_, i), mat; T = Tinit(x[i]))) for i in 1:n]
    G = λ * area / dx
    faces = [connect!(sys, Conductor(ThermalConductance(G)), cells[i], cells[i+1]; name = Symbol(name, :_f, i))
             for i in 1:n-1]
    left_face = _rod_boundary!(sys, left, cells[1], area, λ, dx, Symbol(name, :_left), :left)
    right_face = _rod_boundary!(sys, right, cells[end], area, λ, dx, Symbol(name, :_right), :right)
    all_faces = vcat(left_face === nothing ? Symbol[] : [left_face], faces,
                     right_face === nothing ? Symbol[] : [right_face])
    return (; cells, faces, left_face, right_face, all_faces, x, dx)
end

_rod_boundary!(sys, ::Nothing, cell, area, λ, dx, name, side) = nothing

function _rod_boundary!(sys, res::Symbol, cell, area, λ, dx, name, side)
    return _rod_end!(sys, res, cell, 2λ * area / dx, name, side)
end

function _rod_boundary!(sys, spec::Tuple{Symbol,<:Real}, cell, area, λ, dx, name, side)
    res, h = spec
    h > 0 || throw(ArgumentError("convection coefficient must be positive, got $h"))
    G = 1 / (1 / (h * area) + dx / (2λ * area))
    return _rod_end!(sys, res, cell, G, name, side)
end

function _rod_end!(sys, res, cell, G, name, side)
    a, b = side === :left ? (res, cell) : (cell, res)
    return connect!(sys, Conductor(ThermalConductance(G)), a, b; name)
end
