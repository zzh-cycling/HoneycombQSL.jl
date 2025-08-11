using Random, Test
using QSLQMC

@testset "QSLQMC Tests" begin
    @testset "Greeting Test" begin
        @test QSLQMC.greet() == "Hello World!"
    end

    # Add more tests here as needed
end