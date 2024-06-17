using Feynman
using Test
using Oscar

x=Feynman.sample1(23);

@testset "Feynman.jl" begin
    #sample
    @test Feynman.greet()=="Hello"
    @test Feynman.greet()!="hello"
    
    # Write your tests here.

end
