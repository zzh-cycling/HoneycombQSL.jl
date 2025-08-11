using Test
using LatticeQSL
using Arpack: eigs
using LinearAlgebra

@testset "honeycomb_strings" begin
    m = 3
    n = 4
    xstrings, ystrings, zstrings = honeycomb_strings(m, n)
    
    @test ystrings == [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [10, 11], [12, 13], [14, 15], [16, 17], [18, 19], [20, 21], [22, 23]]
    @test length(ystrings) == m*n
    @test xstrings == [[1, 8], [3, 10], [5, 12], [7, 14], [9, 16], [11, 18], [13, 20], [15, 22], [17, 0], [19, 2], [21, 4], [23, 6]]
    @test length(xstrings) == m*n
    @test zstrings == [[1, 2], [3, 4], [5, 6], [9, 10], [11, 12], [13, 14], [17, 18], [19, 20], [21, 22], [0, 7], [8, 15], [16, 23]]
    @test length(zstrings) == m*n

    m = 2
    n = 4
    xstrings, ystrings, zstrings = honeycomb_strings(m, n, cluster=:hexagonal)

    # @test ystrings == [[1,2],[3,4], [5,6], [7,8], [9,10], [11,12], [13,14], [15,16], [17,18], [19,20], [21,22], [23,24]]
    # @test xstrings == [[2, 9], [4, 11], [6, 13], [8, 15], [10, 17], [12, 19], [14, 21], [16, 23], [18, 1], [20, 3], [22, 5], [24, 7]]
    # @test zstrings == [[2, 3], [4, 5], [6, 7], [10, 11], [12, 13], [14, 15], [18, 19], [20, 21], [22, 23], [1, 8], [9, 16], [17, 24]]
end

@testset "loop_path" begin
    m=3;n=4
    @test loop_path(1,1,m,n) ==[1, 2, 3, 8, 9, 10]
    @test loop_path(1,4,m,n) ==[7, 0, 1, 14, 15, 8]
    @test loop_path(2,4,m,n) ==[15, 8, 9, 22, 23, 16]
    @test loop_path(3,4,m,n) ==[23, 16, 17, 6, 7, 0] 
    @test loop_path(3,1,m,n) ==[17, 18 ,19, 0, 1, 2] 
    @test loop_path(3,2,m,n) ==[19, 20, 21, 2, 3, 4] 
    @test loop_path(3,3,m,n) ==[21, 22, 23, 4, 5, 6] 

    m=3;n=3
    @test loop_path(1,1,m,n) ==[1, 2, 3, 6, 7, 8]
    @test loop_path(1,3,m,n) ==[5, 0, 1, 10, 11, 6]
    @test loop_path(2,3,m,n) ==[11, 6, 7, 16, 17, 12]
    @test loop_path(3,3,m,n) ==[17, 12, 13, 4, 5, 0] 
    @test loop_path(3,1,m,n) ==[13, 14, 15, 0, 1, 2] 
    @test loop_path(3,2,m,n) ==[15, 16, 17, 2, 3, 4]
end

@testset "kitaev_hamiltonian_sparse and loop_op" begin
    m=2;n=2
    H = kitaev_hamiltonian_sparse(m, n)
    @test size(H) == (2^(2*m*n), 2^(2*m*n))
    vals,vecs = eigs(H, nev=5, which=:SR)
    fllis=[loop_op(i, j, m, n) for i in 1:m, j in 1:n]
    gs0=vecs[:,2]
    gs1=vecs[:,3]
    gs2=vecs[:,4]
    loop_val0 = vec([gs0'*fllis[i,j]*gs0 for i in 1:m, j in 1:n])
    loop_val1 = vec([gs1'*fllis[i,j]*gs1 for i in 1:m, j in 1:n])
    loop_val2 = vec([gs2'*fllis[i,j]*gs2 for i in 1:m, j in 1:n])
    @test loop_val0 ≈ ones(m*n)
    @test loop_val1 ≈ ones(m*n)
    @test loop_val2 ≈ ones(m*n)

    # Should all commute
    commutatorlis = [norm(H*fllis[i,j] - fllis[i,j]*H) for i in 1:m, j in 1:n]
    @test all(commutatorlis .< 1e-10)

    # Wilson1, Wilson2, is the reason 4-fold degeneracy, by check the minus sign in "vcat([σx], (-1)^(m-1)*fill(σz, 2*m-2), [σx])", you can find why some cluster do not have 4-fold degeneracy.
    # wilson1, wilson2 = Wilson12(m, n)
    # @test gs0'*wilson1*gs0 ≈ 1
    # @test gs1'*wilson1*gs1 ≈ 1
    # @test gs2'*wilson1*gs2 ≈ 1
    # @test gs0'*wilson2*gs0 ≈ 1
    # @test gs1'*wilson2*gs1 ≈ 1
    # @test gs2'*wilson2*gs2 ≈ 1
end

