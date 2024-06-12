using Feynman
using Test
using Oscar

@testset "Feynman.jl" begin
    #sample
    @test Feynman.greet()=="Hello"
    @test Feynman.greet()!="hello"
    
    # Write your tests here.

end
