function lattice_translations(Lx::Int, Ly::Int)
    """
    For a specified honeycomb lattice (with Lx unit cells in the x-direction, and Ly units cells in the y-direction),
    returns a list of site orderings which are equivalent via translational symmetry (with periodic boundary conditions).

    Parameters
    ----------
    Lx : Int
        Number of unit cells in the x-direction.
    Ly : Int
        Number of unit cells in the y-direction.

    Returns
    -------
    translation_order_list : Vector{Vector{Int}}
        A list of all equivalent site orderings generated via non-zero translations (with periodic boundary conditions).
    """
    unit_cells = zeros(Int, Lx, Ly, 2)
    translation_order_list = Vector{Vector{Int}}()

    for x in 0:(Lx - 1)
        for y in 0:(Ly - 1)
            left_site = (x * (2 * Ly)) + y
            right_site = left_site + Ly
            unit_cells[x + 1, y + 1, :] .= [left_site, right_site]  # Adjust for 1-based indexing
        end
    end

    for j in 0:(Ly - 1)
        for i in 0:(Lx - 1)
            shifted_lattice = circshift(unit_cells, (i, j, 0))

            site_sequence = Int[]

            for x1 in 1:Lx
                for y1 in 1:Ly
                    push!(site_sequence, shifted_lattice[x1, y1, 1])
                end
                for y2 in 1:Ly
                    push!(site_sequence, shifted_lattice[x1, y2, 2])
                end
            end

            push!(translation_order_list, site_sequence)
        end
    end

    return translation_order_list[2:end], unit_cells  # Return all non-zero translations of lattice
end

function mirror_states(state::T, translation_order_list::Vector{Vector{Int}}) where {d, N, T <: DitStr{d, N, Int}}
    """
    For a single state represented by a ternary string, this function finds all equivalent 'mirror' states
    related by symmetry, and their base-10 integer representation. The state with the smallest base-10 integer
    representation is chosen to be the 'representative state'. The number of unique mirror states for a given
    state, n_unique, is also found.

    Parameters
    ----------
    state : String
        A state represented by a ternary string.
    translation_order_list : Vector{Vector{Int}}
        A list of translationally equivalent site orderings.

    Returns
    -------
    sorted_ints : Vector{Int}
        A list of the base-10 integers (in ascending order) representing the equivalent mirror states.
    n_unique : Int
        The number of unique mirror states.
    rep_state : String
        A ternary string corresponding to the state's representative state.
    rep_int : Int
        A base-10 integer corresponding to the state's representative state.
    """
    L = length(state)
    rev_state = reverse(state)

    state_translations = [state]
    state_translation_ints = [state.buf]

    for ordering in translation_order_list
        translated_state = join(rev_state[i] for i in ordering)
        push!(state_translations, reverse(translated_state))
        push!(state_translation_ints, reverse(translated_state).buf)
    end

    state_inversions = []
    state_inversion_ints = []

    for s in state_translations
        push!(state_inversions, reverse(s))
        push!(state_inversion_ints, reverse(s).buf)
    end

    mirror_ints = vcat(state_translation_ints, state_inversion_ints)
    sorted_ints = sort(mirror_ints)
    n_unique = length(unique(sorted_ints))

    rep_int = sorted_ints[1]
    rep_state = T(rep_int)

    return sorted_ints, n_unique, rep_state, rep_int
end

function get_representative_states(L::Int, translation_order_list::Vector{Vector{Int}}, ints_filename::String="kept_ints", states_filename::String="state_map", nunique_filename::String="n_unique_list") where {d, N, T <: DitStr{d,L,Int}}
    """
    For a lattice with L sites and therefore 3^L possible configurations, this function finds all the representative
    states we need to keep, thus reducing the Hilbert space dimension from 3^L to a much smaller number ndim.

    Parameters
    ----------
    L : Int
        The total number of lattice sites.
    translation_order_list : Vector{Vector{Int}}
        A list of translationally equivalent site orderings.
    ints_filename : String
        Filename for the .npy file created which stores the ndim representative states (in base-10).
    states_filename : String
        Filename for the .npy file created which stores the representative states (in base-10) for each of
        the 3^L possible lattice configurations, i.e. a mapping from a state to its representative state.
    nunique_filename : String
        Filename for the .npy file created which stores the number of unique mirror states for each
        representative state.

    Returns
    -------
    kept_ints :: Vector{Int}
        An array of the base-10 integers corresponding to the representative states.
    state_map :: Vector{Int}
        An array containing the representative state (in base-10) for each of the 3^L possible states.
    n_unique_list :: Vector{Int}
        An array containing the number of unique mirror states for each representative state.
    ndim :: Int
        The total number of representative states we keep, i.e. the Hilbert space dimension is reduced
        from 3^L to ndim.
    """
    kept_ints = Int[]
    n_unique_list = Int[]
    state_map = fill(-1, 3^L)
    index = -1

    for i in 0:(3^L - 1)
        if i % 500000 == 0
            @printf("Finding representative for state: %d\n", i)
        end

        if state_map[i + 1] != -1  # Adjust for 1-based indexing
            continue
        end

        index += 1
        push!(kept_ints, i)

        mirror_ints, n_unique, rep_state, rep_int = mirror_states(T(i), translation_order_list)

        for j in mirror_ints
            state_map[j + 1] = index  # Adjust for 1-based indexing
        end

        push!(n_unique_list, n_unique)
    end

    ndim = length(kept_ints)

    kept_ints = collect(kept_ints)
    n_unique_list = collect(n_unique_list)

    # Save the arrays to files
    # You can use your preferred method to save these, e.g., using `JLD` or `Serialization` package
    # Example: JLD.save(ints_filename * ".jld2", "kept_ints", kept_ints)
    # Example: JLD.save(states_filename * ".jld2", "state_map", state_map)
    # Example: JLD.save(nunique_filename * ".jld2", "n_unique_list", n_unique_list)

    return kept_ints, state_map, n_unique_list, ndim
end

function get_neighbors(Lx::Int, Ly::Int, bond::Symbol)
    """
    For a given bond-direction in the honeycomb lattice (x, y, or z), this finds all the pairs of sites
    along this type of bond.

    Parameters
    ----------
    Lx : Int
        Number of unit cells in the x-direction.
    Ly : Int
        Number of unit cells in the y-direction.
    bond : String
        A string which is either "x", "y", or "z" corresponding to the desired honeycomb bond direction.

    Returns
    -------
    neighbors : Vector{Tuple{Int, Int}}
        A list of tuples of site pairs which each constitute a bond of the chosen type.
    """
    # Note: Total number of sites in lattice: L = Lx * Ly * 2
    unit_cells = zeros(Int, Lx, Ly, 2)

    for x in 0:(Lx - 1)
        for y in 0:(Ly - 1)
            left_site = (x * (2 * Ly)) + y
            right_site = left_site + Ly
            unit_cells[x + 1, y + 1, :] .= [left_site, right_site]  # Adjust for 1-based indexing
        end
    end

    neighbors = Vector{Tuple{Int, Int}}()

    if bond == :x
        for i in 0:(Lx - 1)
            for j in 0:(Ly - 1)
                x_bond_left = unit_cells[i + 1, j + 1, 1]
                x_bond_right = unit_cells[i + 1, j + 1, 2]
                x_pair = (min(x_bond_left, x_bond_right), max(x_bond_left, x_bond_right))
                push!(neighbors, x_pair)
            end
        end
    elseif bond == :y
        for i in 0:(Lx - 1)
            for j in 0:(Ly - 1)
                y_bond_left = unit_cells[i + 1, j + 1, 2]
                y_bond_right = unit_cells[(i + 1) % Lx + 1, j + 1, 1]  # Wrap around in x-direction
                y_pair = (min(y_bond_left, y_bond_right), max(y_bond_left, y_bond_right))
                push!(neighbors, y_pair)
            end
        end
    elseif bond == :z
        for i in 0:(Lx - 1)
            for j in 0:(Ly - 1)
                z_bond_left = unit_cells[i + 1, (j + 1) % Ly + 1, 1]  # Wrap around in y-direction
                z_bond_right = unit_cells[i + 1, j + 1, 2]
                z_pair = (min(z_bond_left, z_bond_right), max(z_bond_left, z_bond_right))
                push!(neighbors, z_pair)
            end
        end
    end

    return neighbors
end


# # Example parameters
# Lx = 3
# Ly = 2

# # Get the lattice translations
# translations = lattice_translations(Lx, Ly)

# # Example state in ternary
# state = "120"

# # Get the representative states
# kept_ints, state_map, n_unique_list, ndim = get_representative_states(4, translations)

# # Get neighbors for a specified bond direction
# neighbors_x = get_neighbors(Lx, Ly, "x")
# neighbors_y = get_neighbors(Lx, Ly, "y")
# neighbors_z = get_neighbors(Lx, Ly, "z")

# println("Neighbors in x-direction: ", neighbors_x)
# println("Neighbors in y-direction: ", neighbors_y)
# println("Neighbors in z-direction: ", neighbors_z)
