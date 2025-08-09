using Test
# include("jl")   # Adjust this import based on your module structure
# include("hamiltonian_functions.jl") # Adjust this import based on your module structure
# include("ground_state_functions.jl") # Adjust this import based on your module structure

function test_ground_state()
    L = 12
    Lx = 3
    Ly = 2

    x_neighbors = get_neighbors(Lx, Ly, bond="x")
    y_neighbors = get_neighbors(Lx, Ly, bond="y")
    z_neighbors = get_neighbors(Lx, Ly, bond="z")

    translation_order_list = lattice_translations(Lx, Ly)

    kept_ints, state_map, n_unique_list, ndim = get_representative_states(
        L, translation_order_list,
        ints_filename="kept_ints_L12",
        states_filename="state_map_L12",
        nunique_filename="n_unique_list_L12"
    )

    Hxx = hamiltonian_functions.construct_Hxx(L, x_neighbors, kept_ints, state_map, n_unique_list, filename="Hxx_L12")
    Hyy = hamiltonian_functions.construct_Hyy(L, y_neighbors, kept_ints, state_map, n_unique_list, filename="Hyy_L12")
    Hzz = hamiltonian_functions.construct_Hzz(L, z_neighbors, kept_ints, filename="Hzz_L12")
    HD = hamiltonian_functions.construct_HD(L, kept_ints, state_map, n_unique_list, filename="HD_L12")

    Kx = 1.0
    Ky = 1.0
    Kz = 1.0
    D = 0.0

    H = (Kx * Hxx) + (Ky * Hyy) + (Kz * Hzz) + (D * HD)

    E_0, psi_0 = ground_state_functions.ground_state(H, L)

    # Clean up memory
    empty!(Hxx)
    empty!(Hyy)
    empty!(Hzz)
    empty!(HD)
    empty!(H)
    empty!(psi_0)

    # Collect garbage
    GC.gc()

    @test -0.67 < E_0 < -0.66
end

# Run the ground state test
test_ground_state()