@doc raw"""
    FThe package Feynman generate the Feynman integrand associated to a Feynman diagram (graph) and compute an IBP system using the powerful module-intersection integration-by-parts (IBP) method, suitable for multi-loop and multi-scale Feynman integral reduction. It will provide( soon! ) an application programming interface(API) in OSCAR to use packages NeatIBP, pfd-parallel to make this computation much faster and to solve the reduction problem associated to Feynman integrals completely.The package Feynman is based on the computer algebra system OSCAR and is provided as a package for the Julia programming language.

This package can generate the Feynman integrand associated to a Feynman graph $G$ if $\text{number of internal edges} < \frac{1}{2}l(l+1)+el$.
Here $l$ is the loop number and $e=span<p_1,...,p_{n_{ext}}>$ is the number of linearly independent external momenta of $G$. In the generic case (That we considered in this package), $e=p_{n_{ext}}-1$.

In the case $\text{number of internal edges} < \frac{1}{2}l(l+1)+el$, the package generate the Feynman integrand associated to a larger Feynman graph $G'$ so that $G$ is a subgraph of  $G'$. User can set appropriately the denomiator powers to zero to obtain the Baikov representation of $G$.

"""
module Feynman
using Oscar

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
export printIBP
export sample1
export Net
export simple_graph
export labeledgraph
export computeM1
export IBP


 include("Functions.jl")
 #include("DataTypes.jl")
# Write your package code here.



end
