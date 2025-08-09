using HoneycombQSL, SparseArrays, LinearAlgebra, Arpack


fllis=[flux(i, j, 3, 3) for i in 1:3, j in 1:3]
H = kitaev_hamiltonian_sparse(3, 3)
@time vals,vecs = eigs(H, nev=5, which=:SR)
save("./Kitaev_eigenm3n3.jld", "vals",vals, "vecs", vecs)

gs=vecs[:,1]
gs1=vecs[:,2]
gs2=vecs[:,3]
vec4 = vecs[:,4]
flux_val = [gs'*fllis[i,j]*gs for i in 1:3, j in 1:3]
[gs2'*fllis[i,j]*gs2 for i in 1:3, j in 1:3]
[gs1'*fllis[i,j]*gs1 for i in 1:3, j in 1:3]
[vec4'*fllis[i,j]*vec4 for i in 1:3, j in 1:3]
using JLD
data=load("/Users/cycling/Documents/Julia/CM_model/myKitaevHoneycomb/Kitaev_eigenm3n4.jld", "vals","vecs")

data=load("/Users/cycling/Documents/Julia/CM_model/myKitaevHoneycomb/Kitaev_eigenm3n3.jld", "vals","vecs")

commutatorlis = [norm(H*fllis[i,j] - fllis[i,j]*H) for i in 1:3, j in 1:3]
@test all(commutatorlis .< 1e-10)

Heven = kitaev_hamiltonian_sparse(2, 4)
valseven,vecseven = eigs(Heven, nev=10, which=:SR)
save("./Kitaev_eigenm2n4.jld", "vals",valseven, "vecs", vecseven)
# Nothing that only even even number boundary sites have 4-fold degeneracy. if 3*3, 3*4, 2*2 do not have 4-fold degeneracy.