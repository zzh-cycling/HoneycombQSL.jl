using HoneycombQSL, SparseArrays, LinearAlgebra, Arpack
using JLD

m=3
n=4
H = kitaev_hamiltonian_sparse(m, n)
@time vals,vecs = eigs(H, nev=5, which=:SR)
save("./exm/data/Kitaev_eigenm$(m)n$(n).jld", "vals",vals, "vecs", vecs)
# Nothing that only even even number boundary sites have 4-fold degeneracy. if 3*3, 3*4, 2*2 do not have 4-fold degeneracy.


function check_Wp(m::Int,n::Int)
    vals, vecs = load("./exm/data/Kitaev_eigenm$(m)n$(n).jld", "vals","vecs")
    @show vals
    fllis = [loop_op(i, j, m, n) for i in 1:m, j in 1:n]
    loop_val_lis = [vec([vecs[:,k]'*fllis[i,j]*vecs[:,k] for i in 1:m, j in 1:n]) for k in 1:size(vecs, 2)]

    return loop_val_lis
end

loop_val_lis24 = check_Wp(2,4)
loop_val_lis33 = check_Wp(3,3)