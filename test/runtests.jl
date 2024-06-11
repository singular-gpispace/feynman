using Feynman
using Test
using Oscar

x=Feynman.sample1(23);

#-----Example for printGraph
G1=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
print(G1.edges)
Feynman.printGraph(G1);
G1
#------Example for printLabeledGraph
var=["x","y","z","p","q","r"];
R, (x,y,z,p,q,r)=polynomial_ring(QQ,var);
G2=Feynman.labeledgraph([1,2,3,4],[(1,3),(1,2),(1,2),(2,4),(3,4),(3,4)],R,var,R,[3,4],R,[[y+1,R(2)] [R(3),r+3]]);
Feynman.printLabeledGraph(G2);

#-----Example for labelGraph
G3=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G4=Feynman.labelGraph(G3,0);
Feynman.printLabeledGraph(G4);
G4

#-----Example for balancingIdeal
G3=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G=Feynman.labelGraph(G3,0);
I=Feynman.balancingIdeal(G);
I

#---------

#---------

#------------------------------------------------------------------
#-----Example for eliminateVariables
G3=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G4=Feynman.labelGraph(G3,0);
G5=Feynman.eliminateVariables(G4);
G5.elimvar
G5.labels

#-----Example for removeVariable
R,v=polynomial_ring(QQ,"p"=>(1:3));
Feynman.removeVariable(R,[2])

#-----Example for removeParameter
R,v,w=polynomial_ring(QQ,"p"=>(1:5),"q"=>(1:6));

I=ideal(R,[v[1],v[2],v[3],v[4],v[5],w[1]])
u=complement_of_prime_ideal(I);
S,iso=localization(R,u);
S
Feynman.removeParameter(S,[2])
#---------------
#-----Example for feynmanDenominators
G=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G1=Feynman.labelGraph(G,0);
Gelim=Feynman.eliminateVariables(G1);
I=Feynman.feynmanDenominators(Gelim);
I

#-----Example for propagators
G=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G=Feynman.simple_graph([1,2,3,4,5,6],[(1,2),(3,6),(4,5),(1,6),(2,3),(5,6),(3,4),1,2,5,4]);
G1=Feynman.labelGraph(G,0);
Gelim=Feynman.eliminateVariables(G1);
Gelim=Feynman.removeElimVars(Gelim);
Gelim.labels;
I=Feynman.propagators(Gelim);



#-----------------------------
#-----------------------------

#-----Example for ISP
G=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G=Feynman.simple_graph([1,2,3,4,5,6],[(1,2),(3,6),(4,5),(1,6),(2,3),(5,6),(3,4),1,2,5,4]);
G=Feynman.simple_graph([1,2,3,4,5,6],[(1,7),(1,2),(2,6),(6,3),(3,4),(7,4),(7,5),(5,6),1,2,3,4,5]);
G1=Feynman.labelGraph(G,0);
Gelim=Feynman.eliminateVariables(G1);
Feynman.ISP(Gelim)

#--------------------------------------

#--------------------------------------
#-----Example for removeElimVars
G=Feynman.simple_graph([1,2,3,4,5,6],[(1,2),(3,6),(4,5),(1,6),(2,3),(5,6),(3,4),1,2,5,4]);
G=Feynman.simple_graph([1,2,3,4,5,6],[(1,7),(1,2),(2,6),(6,3),(3,4),(7,4),(7,5),(5,6),1,2,3,4,5]);
G1=Feynman.labelGraph(G,0);
Gelim=Feynman.eliminateVariables(G1);
Gelim.labels
G=Feynman.removeElimVars(Gelim);
G
#-----------------------------------
#-----------------------------------
#-----Example for computeBaikovMatrix
#Example 1
G=Feynman.simple_graph([1,2,3,4,5,6],[(1,2),(3,6),(4,5),(1,6),(2,3),(5,6),(3,4),1,2,5,4]);
G1=Feynman.labelGraph(G,0);
Gelim=Feynman.eliminateVariables(G1);
propagators(Gelim)
Gn=Feynman.removeElimVars(Gelim);
Gn.labels
propagators(Gn)
G=Feynman.computeBaikovMatrix(Gn);
G.baikovmatrix

#Example 2
G=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G1=Feynman.labelGraph(G,0);
Gelim=Feynman.eliminateVariables(G1);
Gn=Feynman.removeElimVars(Gelim);
G=Feynman.computeBaikovMatrix(Gn);
G.baikovmatrix
printLabeledGraph(G)

#Example 3
G=Feynman.simple_graph([1,2,3,4,5,6],[(1,2),(1,3),(2,4),(4,5),(5,6),(2,6),(3,6),1,2,3,4]);
G=Feynman.labelGraph(G,0);
G=Feynman.eliminateVariables(G);
G=removeElimVars(G);
G=computeBaikovMatrix(G);
G.baikovmatrix
#Example 4
G=Feynman.simple_graph([1,2,3,4,5,6],[(1,2),(2,4),(3,1),(4,3),(6,4),(3,5),(5,6),1,2,6,5]);
G=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G=Feynman.labelGraph(G,0);
G=Feynman.eliminateVariables(G);
G=Feynman.removeElimVars(G);
G=Feynman.computeBaikovMatrix(G);
G.baikovmatrix
#-----------
G=Feynman.simple_graph([1,2,3,4,5,6,7],[(1,7),(1,2),(2,6),(6,3),(3,4),(7,4),(7,5),(5,6),1,2,3,4,5]);
G=Feynman.labelGraph(G,0);
G=Feynman.eliminateVariables(G);
G=Feynman.removeElimVars(G);
G=Feynman.computeBaikovMatrix(G);
G.baikovmatrix

#-------
#-------
@testset "Feynman.jl" begin
    #sample
    @test Feynman.greet()=="Hello"
    @test Feynman.greet()!="hello"
    @test x.a==23
    # Write your tests here.

end
