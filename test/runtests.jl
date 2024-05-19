using Feynman
using Test
using Oscar

x=Feynman.sample1(23);

#-----Example for printGraph
G1=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
print(G1.edges)
Feynman.printGraph(G1);

#------Example for printLabeledGraph
var=["x","y","z","p","q","r"];
R, (x,y,z,p,q,r)=polynomial_ring(QQ,var);
invlex(R)
x+y
G2=Feynman.labeledgraph([1,2,3,4],[(1,3),(1,2),(1,2),(2,4),(3,4),(3,4)],R,var,R,[3,4],R,[[y+1,R(2)] [R(3),r+3]]);
Feynman.printLabeledGraph(G2);

#-----Example for labelGraph
G3=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G4=Feynman.labelGraph(G3,0);
Feynman.printLabeledGraph(G4);

#-----Example for balancingIdeal
G3=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G4=Feynman.labelGraph(G3,0);
I=Feynman.balancingIdeal(G4);
I

#-----Example for eliminateVariables
G3=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G4=Feynman.labelGraph(G3,0);
G5=Feynman.eliminateVariables(G4);
G5

#-----Example for removeVariable
R,v=polynomial_ring(QQ,"p"=>(1:3));
Feynman.removeVariable(R,2)

#-----Example for removeParameter
R,v=polynomial_ring(QQ,"p"=>(1:4));
I=ideal(R,[v[1],v[2]]);
u=complement_of_prime_ideal(I);
S,iso=localization(R,u);
Feynman.removeParameter(S,2)

#-----Example for feynmanDenominators
G=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G1=Feynman.labelGraph(G,0);
Gelim=Feynman.eliminateVariables(G1);
I=Feynman.feynmanDenominators(Gelim);
I

#-----Example for propagators
G=Feynman.simple_graph([1,2,3,4],[(1,3),(1,2),(2,4),(3,4),1,2,3,4]);
G1=Feynman.labelGraph(G,0);
Gelim=Feynman.eliminateVariables(G1);
I=Feynman.propagators(Gelim);
I
#--------
G=Gelim;
S=G.over;
RP=G.overpoly;
J=Feynman.propagators(G);
el=G.elimvar;

v=[];
for i in 1:length(G.edges)
    if length(G.edges[i])==1
        push!(v,G.labels[i]);
    end   
end
infedges=ideal(S,v);

T=hom(S,RP,gens(RP));

v=[];
for i in 1:ngens(J)
push!(v,T(J[i])); 
end
J=ideal(RP,v);

v=[];
for i in 1:ngens(infedges)
push!(v,T(infedges[i])); 
end
infedges=ideal(RP,v);

J=J+infedges^2;
J
v=[];
for i in 1:length(el)
push!(v,T(el[i])); 
end

el=v;

for i in 1:length(el)
J=J+ideal(RP,el[i]);    
end

J=standard_basis(J);
J
L=kbase(J,2);
ideal(S,L);

#--------
@testset "Feynman.jl" begin
    #sample
    @test Feynman.greet()=="Hello"
    @test Feynman.greet()!="hello"
    @test x.a==23
    # Write your tests here.

end
