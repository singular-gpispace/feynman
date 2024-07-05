@doc raw"""
    Feynman is a package for computing integration-by-part identities (IBPs) of a Feynman Integral associated to Feynman graph using module intesecation method. 
    This package also provides an interface of Oscar to use the packages NeatIBP developed using Singular and GPI-Space. 

"""
module Feynman
using Oscar

export greet
export printNet
export printGraph
export printLabeledGraph
export labelGraph
export balancingIdeal
export substituteGraph
export eliminateVariables
export removeVariable
export removeParameter
export feynmanDenominators
export propagators
export ISP
export removeElimVars
export computeBaikovMatrix
export makePoly
export removeVariableLocal
export computeIBP

export sample1
export Net
export simple_graph
export labeledgraph
 include("Functions.jl")
 #include("DataTypes.jl")
# Write your package code here.



end
