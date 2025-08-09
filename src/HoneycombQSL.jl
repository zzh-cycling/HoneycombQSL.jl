module HoneycombQSL

using LinearAlgebra, SparseArrays, Random, Arpack, BitBasis
export flux_path, kitaev_hamiltonian_sparse, Wilson12, flux

include("EDdemo.jl")
end
