module HoneycombQSL

using LinearAlgebra, SparseArrays, Random, Arpack, BitBasis
export flux_path, kitaev_hamiltonian_sparse, Wilson12, flux, honeycomb_strings

include("EDdemo.jl")
end
