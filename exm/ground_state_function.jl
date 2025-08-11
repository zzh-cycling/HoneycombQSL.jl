using LatticeQSL
using SparseArrays, LinearAlgebra, Arpack

# Define the Hamiltonian H, total sites L, etc.
H = sprand(100, 100, 0.1)  # Example sparse matrix
L = 100

# Calculate ground state and energy
E_0, psi_0 = ground_state(H, L)

# Assuming HD is defined
HD = sprand(100, 100, 0.1)  # Example sparse matrix for single-ion anisotropy

# Calculate <(Sx + Sy + Sz)^2>
Ss = s_squared(HD, L, psi_0)

# For entanglement entropy, assume state_map and n_unique_list are defined
state_map = rand(1:3^L, 3^L)
n_unique_list = rand(1:10, 3^L)

entropy = entanglement_entropy(L, psi_0, state_map, n_unique_list)

println("Entanglement entropy: ", entropy)
