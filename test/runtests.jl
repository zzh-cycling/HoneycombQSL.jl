using LatticeQSL
using Test

@testset "EDdemo.jl" begin
    # Write your tests here.
    include("EDdemo.jl")
end

@testset "Basis.jl" begin
    # Write your tests here.
    include("Basis.jl")
end

@testset "symmetry" begin
    include("test_symmetry.jl")
end

@testset "test_ground_state.jl" begin
    include("test_ground_state.jl")
end