# Model graph (DESIGN §4.1, D-17): a bipartite MetaGraph whose vertices are nodes (bodies,
# reservoirs) and elements (conductors, streams, pumps, flow boundaries, outlets), and whose
# edges are ports carrying the port role.
# The graph is for building, validation and topology queries; `compile` lowers it to flat,
# type-stable data and the graph is never touched in the right-hand side.

"A port: the edge between an element and a node. `role` is :a, :b, or :h (a conductor's third thermal node)."
struct Port
    role::Symbol
end

const ModelGraph = typeof(MetaGraph(Graph(); label_type = Symbol, vertex_data_type = Any, edge_data_type = Port))

"""
    ThermoSystem()

An empty model. Add bodies and reservoirs with `add!`, elements with `connect!`, then `compile`.
"""
mutable struct ThermoSystem
    graph::ModelGraph
    nodes::Vector{Symbol}       # insertion order
    elements::Vector{Symbol}    # insertion order
end

ThermoSystem() = ThermoSystem(MetaGraph(Graph(); label_type = Symbol, vertex_data_type = Any, edge_data_type = Port),
                              Symbol[], Symbol[])

isnode(x) = x isa Body || x isa Reservoir || x isa Junction

function _check_new_label(sys::ThermoSystem, name::Symbol)
    haskey(sys.graph, name) && throw(ArgumentError("a vertex named :$name already exists"))
    return nothing
end

"""
    add!(sys, node) -> Symbol

Add a node (a `Body`, `Reservoir` or `Junction`); returns its label.
"""
function add!(sys::ThermoSystem, x::Union{Body,Reservoir,Junction})
    _check_new_label(sys, x.name)
    sys.graph[x.name] = x
    push!(sys.nodes, x.name)
    return x.name
end

function _check_node(sys::ThermoSystem, label::Symbol)
    haskey(sys.graph, label) || throw(ArgumentError("no vertex named :$label"))
    isnode(sys.graph[label]) || throw(ArgumentError(":$label is an element, not a body or reservoir"))
    return nothing
end

"""
    connect!(sys, conductor, a, b; name) -> Symbol

Join nodes `a` and `b` with a conductor. Positive flux runs from `a` to `b`. Parallel
conductors between the same pair are allowed. If the conductor's dissipation policy names
a third thermal node (`Into(h)`, `InternalNode(h = …)`), that node gets a port `:h`.
"""
function connect!(sys::ThermoSystem, el::Conductor, a::Symbol, b::Symbol;
                  name::Symbol = Symbol(:conductor_, length(sys.elements) + 1))
    _check_node(sys, a)
    _check_node(sys, b)
    a == b && throw(ArgumentError("a conductor needs two distinct nodes, got :$a twice"))
    h = third_node(el.dissipation)
    if h !== nothing
        _check_node(sys, h)
        h in (a, b) && throw(ArgumentError("the third node must differ from the conductor's ends; " *
                                           "use Split(w) to put the heat into an end"))
    end
    _check_new_label(sys, name)
    sys.graph[name] = el
    sys.graph[name, a] = Port(:a)
    sys.graph[name, b] = Port(:b)
    h === nothing || (sys.graph[name, h] = Port(:h))
    push!(sys.elements, name)
    return name
end

_default_name(sys::ThermoSystem, el) = Symbol(lowercase(string(nameof(typeof(el)))), :_, length(sys.elements) + 1)

function _add_element!(sys::ThermoSystem, el, name::Symbol, ports::Pair{Symbol,Symbol}...)
    _check_new_label(sys, name)
    sys.graph[name] = el
    for (node, role) in ports
        sys.graph[name, node] = Port(role)
    end
    push!(sys.elements, name)
    return name
end

"""
    connect!(sys, stream_or_pump, a, b; name) -> Symbol

Join nodes `a` and `b` with a `Stream` (positive ṁ from a to b) or an `IdealPump` (from a to b).
"""
function connect!(sys::ThermoSystem, el::Union{Stream,IdealPump}, a::Symbol, b::Symbol;
                  name::Symbol = _default_name(sys, el))
    _check_node(sys, a)
    _check_node(sys, b)
    a == b && throw(ArgumentError("a $(nameof(typeof(el))) needs two distinct nodes, got :$a twice"))
    return _add_element!(sys, el, name, a => :a, b => :b)
end

"""
    connect!(sys, boundary_or_outlet, node; name) -> Symbol

Attach a `FlowBoundary` (inflow into `node`) or a `PassThroughOutlet` (outflow from `node`).
"""
function connect!(sys::ThermoSystem, el::Union{FlowBoundary,PassThroughOutlet}, node::Symbol;
                  name::Symbol = _default_name(sys, el))
    _check_node(sys, node)
    return _add_element!(sys, el, name, node => :a)
end

"""
    connect!(sys, source, body; name) -> Symbol

Attach a `Source`, `GravityField`, `HeatSource`, `NuclearSource` or `Reaction` to a body.
"""
function connect!(sys::ThermoSystem, el::Union{Source,GravityField,Reaction,NuclearSource,HeatSource}, node::Symbol;
                  name::Symbol = _default_name(sys, el))
    _check_node(sys, node)
    return _add_element!(sys, el, name, node => :a)
end

"""
    connect!(sys, converter, nodes...; name) -> Symbol

Join a `Converter`'s ports, in order, to distinct nodes (port k is role `:p<k>`).
"""
function connect!(sys::ThermoSystem, el::Converter{N}, nodes::Symbol...;
                  name::Symbol = _default_name(sys, el)) where {N}
    length(nodes) == N || throw(ArgumentError("this converter has $N ports, got $(length(nodes)) nodes"))
    allunique(nodes) || throw(ArgumentError("converter ports must be on distinct nodes"))
    foreach(n -> _check_node(sys, n), nodes)
    return _add_element!(sys, el, name, (n => Symbol(:p, k) for (k, n) in enumerate(nodes))...)
end

"""
    connect!(sys, transformer, a, b; name) -> Symbol

Join a `Transformer` between nodes `a` and `b` (P_a = r P_b).
"""
function connect!(sys::ThermoSystem, el::Transformer, a::Symbol, b::Symbol; name::Symbol = _default_name(sys, el))
    _check_node(sys, a)
    _check_node(sys, b)
    a == b && throw(ArgumentError("a transformer needs two distinct nodes"))
    return _add_element!(sys, el, name, a => :a, b => :b)
end

"""
    connect!(sys, switched, nodes...; name) -> Symbol

Connect a `Switched` element exactly as its inner element would be connected.
"""
function connect!(sys::ThermoSystem, el::Switched, nodes::Symbol...; name::Symbol = _default_name(sys, el.element))
    label = connect!(sys, el.element, nodes...; name)
    sys.graph[label] = el
    return label
end

"Label of the node on port `role` of element `el`."
function port_node(sys::ThermoSystem, el::Symbol, role::Symbol)
    for nb in neighbor_labels(sys.graph, el)
        sys.graph[el, nb].role == role && return nb
    end
    throw(ArgumentError("element :$el has no port $role"))
end
