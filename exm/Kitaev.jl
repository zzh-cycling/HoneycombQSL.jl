using HoneycombQSL, SparseArrays, LinearAlgebra, Arpack
using JLD

H = kitaev_hamiltonian_sparse(3, 3)
@time vals,vecs = eigs(H, nev=20, which=:SR)
save("./Kitaev_eigenm3n3.jld", "vals",vals, "vecs", vecs)



vals24, vecs24 = load("./exm/data/Kitaev_eigenm2n4.jld", "vals","vecs")
gs0 = vecs24[:,1]
gs1 = vecs24[:,2]
gs2 = vecs24[:,3]
gs3 = vecs24[:,4]
fllis = [loop_op(i, j, 2, 4) for i in 1:2, j in 1:4]
loop_val0 = vec([gs0'*fllis[i,j]*gs0 for i in 1:2, j in 1:4])
loop_val1 = vec([gs1'*fllis[i,j]*gs1 for i in 1:2, j in 1:4])
loop_val2 = vec([gs2'*fllis[i,j]*gs2 for i in 1:2, j in 1:4])
loop_val3 = vec([gs3'*fllis[i,j]*gs3 for i in 1:2, j in 1:4])


vals33, vecs33 = load("./exm/data/Kitaev_eigenm3n3.jld", "vals","vecs")
gs0 = vecs33[:,1]
gs1 = vecs33[:,2]
gs2 = vecs33[:,3]
gs3 = vecs33[:,4]
fllis = [loop_op(i, j, 3, 3) for i in 1:3, j in 1:3]
loop_val0 = vec([gs0'*fllis[i,j]*gs0 for i in 1:3, j in 1:3])
loop_val1 = vec([gs1'*fllis[i,j]*gs1 for i in 1:3, j in 1:3])
loop_val2 = vec([gs2'*fllis[i,j]*gs2 for i in 1:3, j in 1:3])
loop_val3 = vec([gs3'*fllis[i,j]*gs3 for i in 1:3, j in 1:3])


vals34, vecs34 = load("./exm/data/Kitaev_eigenm3n4.jld", "vals","vecs")
gs0 = vecs34[:,1]
gs1 = vecs34[:,2]
gs2 = vecs34[:,3]
gs3 = vecs34[:,4]
gs4 = vecs34[:,5]
fllis = [loop_op(i, j, 3, 4) for i in 1:3, j in 1:4]
# loop_val0, loop_val1 all 1. But others fails.
loop_val0 = vec([gs0'*fllis[i,j]*gs0 for i in 1:3, j in 1:4])
loop_val1 = vec([gs1'*fllis[i,j]*gs1 for i in 1:3, j in 1:4])
loop_val2 = vec([gs2'*fllis[i,j]*gs2 for i in 1:3, j in 1:4])
loop_val3 = vec([gs3'*fllis[i,j]*gs3 for i in 1:3, j in 1:4])
loop_val4 = vec([gs4'*fllis[i,j]*gs4 for i in 1:3, j in 1:4])   

Heven = kitaev_hamiltonian_sparse(3, 3)
valseven,vecseven = eigs(Heven, nev=10, which=:SR)
save("./exm/data/Kitaev_eigenm3n3.jld", "vals",valseven, "vecs", vecseven)
# Nothing that only even even number boundary sites have 4-fold degeneracy. if 3*3, 3*4, 2*2 do not have 4-fold degeneracy.