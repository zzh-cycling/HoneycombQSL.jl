module HoneycombQSL

using LinearAlgebra, SparseArrays, Random, Arpack, BitBasis

export kitaev_hamiltonian_sparse, wilson12, loop_path, loop_op, honeycomb_strings
export honeycomb_basis, kitaev_honeycomb_ham
export kagome_strings, kagome_basis, Heisenberg_hamiltonian_sparse, loop_map

include("EDdemo.jl")
include("Basis.jl")
end
