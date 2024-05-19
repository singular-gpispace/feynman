using Oscar
export sample1, Net, simple_graph, labeledgraph

mutable struct sample1
    a::Int
end
mutable struct simple_graph
    vertices :: Vector{Int64}
    edges :: Vector
end

mutable struct Net
    rows :: Vector
end



mutable struct labeledgraph
    vertices :: Vector{Int64}
    edges :: Vector
    over :: Ring
    labels :: Vector
    overpoly ::Ring
    elimvar ::Vector
    baikovover :: Ring
    baikovmatrix :: Matrix{RingElem}
end