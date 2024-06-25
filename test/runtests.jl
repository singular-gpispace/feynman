using Feynman
using Test
using Oscar

x=Feynman.sample1(23);

#=Working place
G=Feynman.simple_graph([1,2,3,4,5,6,7],[(6,1),(6,4),(1,2),(3,7),(4,3),(2,7),(5,6),(7,5),1,2,3,4,5]);
G=Feynman.simple_graph([1,2,3,4,5,6],[(6,1),(4,6),(1,2),(3,5),(4,3),(2,5),(5,6),1,2,3,4]);
G=Feynman.simple_graph([1,2,3,4],[(1,2),(2,3),(3,4),(4,1),1,2,3,4]);
G=Feynman.labelGraph(G,0);
printLabeledGraph(G)
G=Feynman.eliminateVariables(G);
G=Feynman.removeElimVars(G);
G=Feynman.computeBaikovMatrix(G);
B=G.baikovmatrix

for i in 1:6 
    print("x(",i,"5)=",B[i,5]," ,");
end
for i in 1:6 
    print("x(5,",i,")=",B[5,i]," ,");
end
for i in 1:5 
    print(" x(",i,",4)=",B[i,4]," ,");
end
for i in 1:5 
    print("x(",i,",5)=",B[i,5]," ,");
end
for i in 1:5 
    print(" x(1,",i,")=",B[1,i]," ,");
end

for i in 1:5 
    print(" x(4,",i,")=",B[4,i]," ,");
end

for i in 1:4 
    for j in 1:4 
        print("(",i,j,")=",B[i,j]," ");
    end
end




#-----------------------------------------------------
#-----------------------------------------------------
=#
#-----------------------------------------------------
@testset "Feynman.jl" begin
    #sample
    @test Feynman.greet()=="Hello"
    @test Feynman.greet()!="hello"
    
    # Write your tests here.

end
