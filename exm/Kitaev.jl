using LatticeQSL, SparseArrays, LinearAlgebra, Arpack
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

function check_excitation(m::Int,n::Int, idx::Vector{Int})
    vals, vecs = load("./exm/data/Kitaev_eigenm$(m)n$(n).jld", "vals","vecs")
    @show vals
    fllis = [loop_op(i, j, m, n) for i in 1:m, j in 1:n]
    loop_val_mat = [vecs[:,idx]'*fllis[i,j]*vecs[:,idx] for i in 1:m, j in 1:n]

    return loop_val_mat
end

# eigvals of 3*3, 2*4, 3*4
# -14.2915  -12.9443 -19.09178532237227
#  -14.2915  -12.9443 -18.98699905558389
#  -14.2915  -12.9443 -18.82203651360451
#  -14.0746  -12.9443 -18.822036513604388
#  -14.0746  -12.7929 -18.822036513604083
#  -14.0746  -12.7929
#  -14.0746  -12.7929
#  -14.0746  -12.7929
#  -14.021   -12.5851
#  -14.021   -12.5851

loop_val_lis24 = check_Wp(2,4)
loop_val_lis33 = check_Wp(3,3)

loop_val_mat = check_excitation(2,4, collect(5:8))
