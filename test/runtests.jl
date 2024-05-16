using Feynman
using Test
using Oscar
x=Feynman.sample1(23);
@testset "Feynman.jl" begin
    @test Feynman.greet()=="Hello"
    @test Feynman.greet()!="hello"
    @test x.a==23
    # Write your tests here.
end
