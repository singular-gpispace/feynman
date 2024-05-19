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

 include("Functions.jl")
 include("DataTypes.jl")
# Write your package code here.



end
