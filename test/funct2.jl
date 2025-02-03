@testset "funct2.jl" begin

    G=labelGraph(simple_graph([1,2],[(1,2),(2,1),1,2]),0);
    G1=eliminateVariables(G);
    P=G1.overpoly;
    S=G1.over;
    p=gens(P);
    q=[p[3],p[4]];
    I1=ideal([q[1]^2,(p[1]+q[1])^2]);
    I2=ideal([q[1]^2,q[1]^2+2*q[1]*p[1]]);
    H=hom(S,P,gens(P));
    u=[];
    
  for i in 1:2
    push!(u,H(I1[i])); 
    
  end
    I3=ideal(S,u);
    
   
    @testset "feynmanDenominators test" begin
        @test   I3==feynmanDenominators(G1)
    end
    Gelim=eliminateVariables(G);
    G2=removeElimVars(Gelim);

    P=G2.overpoly;
    S=G2.over;
    p=gens(P);
    q=[p[3]];
    I2=ideal([q[1]^2,q[1]^2+2*q[1]*p[1]]);
    H=hom(S,P,gens(P));
    u=[];
    
    for i in 1:2
        push!(u,H(I2[i])); 
    end

    I4=ideal(S,u);

    @testset "propagators test" begin
        @test   I4==propagators(G2)
    end

    G5=labelGraph(simple_graph([1,2],[(1,2),(2,1),1,2]),0);
    G6=eliminateVariables(G5);
    I5=ideal(G6.over,[0]);
    
    @testset "ISP test" begin
        @test   I5==ISP(G6)
    end

    G7=labelGraph(simple_graph([1,2],[(1,2),(2,1),1,2]),0);
    G8=eliminateVariables(G5);
    p=gens(G6.overpoly);
    q=[p[3],p[4]];
    
    u=[q[1],p[1]+q[1],p[1],-p[1]];
    v=[];
    G8=removeElimVars(G8);
    @testset "removeElimVars test" begin
        @test  u==G8.labels
        @test   v==G8.elimvar
    end

    G9=simple_graph([1,2,3,4,5,6],[(6,1),(4,6),(1,2),(3,5),(4,3),(2,5),(5,6),1,2,3,4]);
    G10=labelGraph(G9,0);
    G11=eliminateVariables(G10);
    G12=removeElimVars(G11);
    G13=computeBaikovMatrix(G12);
    G14=computeBaikovMatrix(G9);
    M1=G13.baikovmatrix;
    M2=G14.baikovmatrix;
    z1=gens(G13.baikovover);
    z2=gens(G14.baikovover)
    @testset "computeBaikovMatrix test" begin
        @test z1[5]-z1[6]==2*M1[4,1]
        @test z2[5]-z2[6]==2*M2[4,1]
    end



end