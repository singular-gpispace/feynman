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
 include("Functions.jl")
 include("DataTypes.jl")
# Write your package code here.



end
