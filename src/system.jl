# Model graph (DESIGN §4.1, D-17): a bipartite MetaGraph whose vertices are nodes (bodies,
# reservoirs) and elements (conductors), and whose edges are ports carrying the port role.
# The graph is for building, validation and topology queries; `compile` lowers it to flat,
# type-stable data and the graph is never touched in the right-hand side.

"A port: the edge between an element and a node. `role` is :a, :b, or :h (third thermal node) for a conductor."
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

isnode(x) = x isa Body || x isa Reservoir

function _check_new_label(sys::ThermoSystem, name::Symbol)
    haskey(sys.graph, name) && throw(ArgumentError("a vertex named :$name already exists"))
    return nothing
end

"""
    add!(sys, body_or_reservoir) -> Symbol

Add a node; returns its label.
"""
function add!(sys::ThermoSystem, x::Union{Body,Reservoir})
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

"Label of the node on port `role` of element `el`."
function port_node(sys::ThermoSystem, el::Symbol, role::Symbol)
    for nb in neighbor_labels(sys.graph, el)
        sys.graph[el, nb].role == role && return nb
    end
    throw(ArgumentError("element :$el has no port $role"))
end
