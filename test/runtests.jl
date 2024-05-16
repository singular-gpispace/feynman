using Feynman
using Test

@testset "Feynman.jl" begin
    @test Feynman.greet()=="Hello"
    @test Feynman.greet()=="hello"
    # Write your tests here.
end
