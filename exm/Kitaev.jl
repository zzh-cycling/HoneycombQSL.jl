using HoneycombQSL, SparseArrays, LinearAlgebra, Arpack
using JLD

H = kitaev_hamiltonian_sparse(3, 3)
@time vals,vecs = eigs(H, nev=5, which=:SR)
save("./Kitaev_eigenm3n3.jld", "vals",vals, "vecs", vecs)

data=load("/Users/cycling/Documents/Julia/CM_model/myKitaevHoneycomb/Kitaev_eigenm3n4.jld", "vals","vecs")

data=load("/Users/cycling/Documents/Julia/CM_model/myKitaevHoneycomb/Kitaev_eigenm3n3.jld", "vals","vecs")


Heven = kitaev_hamiltonian_sparse(2, 4)
valseven,vecseven = eigs(Heven, nev=10, which=:SR)
save("./Kitaev_eigenm2n4.jld", "vals",valseven, "vecs", vecseven)
# Nothing that only even even number boundary sites have 4-fold degeneracy. if 3*3, 3*4, 2*2 do not have 4-fold degeneracy.