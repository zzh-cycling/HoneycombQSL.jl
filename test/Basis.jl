using Test
using HoneycombQSL
using Arpack: eigs
using LinearAlgebra
using BitBasis

@testset "loop map" begin
    m = 3
    n = 4
    T = BitStr{2*m*n}
    state = T(0)
    
    @test loop_map(1, 1, m, n, state) == (T(bit"011000000110000000000000"),  -1)
    @test loop_map(1, 4, m, n, state) == (T(bit"100000011000000100000000"),  -1)
    @test loop_map(2, 4, m, n, state) == (T(bit"000000001000000110000001"),  -1)
    @test loop_map(3, 4, m, n, state) == (T(bit"100000010000000010000001"),  -1)
end


m = 2
n = 3
T = BitStr{2*m*n}
loop_map.(1, 1, m, n, honeycomb_basis(m,n))