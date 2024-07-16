using Feynman
using Test
using Oscar


include("funct.jl")


#------------------------------------------
#------------------------------------------
G=Feynman.simple_graph([1,2,3,4],[(1,4),(1,2),(2,3),(3,4),1,2,3,4]);
G=simple_graph([1,2,3,4,5,6],[(6,1),(4,6),(1,2),(3,5),(4,3),(2,5),(5,6),1,2,3,4]);
G=simple_graph([1,2,3,4,5,6,7],[(6,1),(6,4),(1,2),(3,7),(4,3),(2,7),(5,6),(7,5),1,2,3,4,5]);

G=Feynman.labelGraph(G,0);
G=Feynman.eliminateVariables(G);
G=Feynman.removeElimVars(G);
G=Feynman.computeBaikovMatrix(G);
computeIBP(G,7,5);
