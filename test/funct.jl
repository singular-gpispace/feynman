@testset "funct.jl" begin
    G=simple_graph([],[]);
    G1=simple_graph([1,2],[(1,2)]);

    @testset "labelGraph test" begin
        @test_throws ErrorException labelGraph(G,1)
        @test_throws ErrorException labelGraph(G1,-1)
    end

    
    G2=labelGraph(simple_graph([1,2],[(1,2),(2,1),1,2]),0);
    P=G2.overpoly;
    p=gens(P);
    q=[p[3],p[4]];
    Q=complement_of_prime_ideal(ideal(P,p));
    R,iso=localization(P,Q);
    I=ideal(G2.over,[p[1]+p[2],p[1]+q[1]-q[2],p[2]+q[2]-q[1]]);

    @testset "balancingIdeal test" begin
    
        @test balancingIdeal(G2)==I
    end
    G2=labelGraph(simple_graph([1,2],[(1,2),(2,1),1,2]),0);
    G3=G2;
    
    G3.labels[1]=-p[2];
    @testset "substituteGraph test" begin
       @test substituteGraph(G2,p[1],-p[2])==G3
    end
    
    G2=labelGraph(simple_graph([1,2],[(1,2),(2,1),1,2]),0);
    G4=G2;
    G4.labels=[];
    G4.labels=[q[1],p[1]+q[1],p[1],-p[1]];
    G4.elimvar=[p[2],q[2]];
    
    @testset "eliminateVariables test" begin
    @test G4==eliminateVariables(G2)    
    end
    
    R,p=polynomial_ring(QQ,"p"=>1:4)
    S,p=polynomial_ring(QQ,["p[1]","p[4]"])

    @testset "removeVariable test" begin
        @test_throws ErrorException removeVariable(R,[5])
         @test   S==removeVariable(R,[2,3])
    end
    
    P,p,q=polynomial_ring(QQ,"p"=>1:4,"q"=>1:3)
            Q=complement_of_prime_ideal(ideal(P,p));
            R,iso=localization(P,Q);
            
            S,p,q=polynomial_ring(QQ,"p"=>1:1,"q"=>1:3)
            
            v=[p[1]];
            T=complement_of_prime_ideal(ideal(S,p));
            U,iso=localization(S,T);
    
    @testset "removeParameter test" begin
        @test_throws ErrorException removeVariable(R,[8])
    end

end