module HoneycombQSL

using LinearAlgebra, SparseArrays, Random, Arpack, BitBasis
export kitaev_hamiltonian_sparse, wilson12, loop_path, loop_op, honeycomb_strings

include("EDdemo.jl")
end
