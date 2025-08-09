# Assuming the necessary functions are defined in the appropriate module
# Or you can include the implementations of `ternary`, `ternary_pad`, and `tern_to_base10` here

function SmSm(state::String, j::Int, k::Int)
    """
    Applies spin lowering operator S- to sites j and k, i.e. acts S-(j) S-(k) on a state.

    Parameters
    ----------
    state : String
        A string representing a state in ternary.
    j : Int
        An integer in the range [0, L-1] (where L is the total number of lattice sites) denoting a site j to lower.
    k : Int
        An integer in the range [0, L-1] (where L is the total number of lattice sites) denoting a site k to lower.

    Returns
    -------
    new_state : String
        The ternary representation of the new state after the spin lowering operator S- has acted on sites j and k.
    """
    
    state_rev = reverse(state)
    list_sites = collect(state_rev)

    list_sites[j + 1] = string(parse(Int, list_sites[j + 1]) + 1)  # Adjusting for 1-based indexing
    list_sites[k + 1] = string(parse(Int, list_sites[k + 1]) + 1)

    new_state = join(reverse(list_sites))

    return new_state
end

function SmSp(state::String, j::Int, k::Int)
    """
    Applies spin lowering operator S- to site j, and spin raising operator S+ to site k, i.e. acts S-(j) S+(k) on a state.

    Parameters
    ----------
    state : String
        A string representing a state in ternary.
    j : Int
        An integer in the range [0, L-1] (where L is the total number of lattice sites) denoting a site j to lower.
    k : Int
        An integer in the range [0, L-1] (where L is the total number of lattice sites) denoting a site k to raise.

    Returns
    -------
    new_state : String
        The ternary representation of the new state after the spin lowering operator S- has acted on sites j, and the spin raising operator S+ has acted on site k.
    """
    
    state_rev = reverse(state)
    list_sites = collect(state_rev)

    list_sites[j + 1] = string(parse(Int, list_sites[j + 1]) + 1)
    list_sites[k + 1] = string(parse(Int, list_sites[k + 1]) - 1)

    new_state = join(reverse(list_sites))

    return new_state
end

function SpSm(state::String, j::Int, k::Int)
    """
    Applies spin raising operator S+ to site j, and spin lowering operator S- to site k, i.e. acts S+(j) S-(k) on a state.

    Parameters
    ----------
    state : String
        A string representing a state in ternary.
    j : Int
        An integer in the range [0, L-1] (where L is the total number of lattice sites) denoting a site j to raise.
    k : Int
        An integer in the range [0, L-1] (where L is the total number of lattice sites) denoting a site k to lower.

    Returns
    -------
    new_state : String
        The ternary representation of the new state after the spin raising operator S+ has acted on site j, and the spin lowering operator S- has acted on site k.
    """
    
    state_rev = reverse(state)
    list_sites = collect(state_rev)

    list_sites[j + 1] = string(parse(Int, list_sites[j + 1]) - 1)
    list_sites[k + 1] = string(parse(Int, list_sites[k + 1]) + 1)

    new_state = join(reverse(list_sites))

    return new_state
end

function SpSp(state::String, j::Int, k::Int)
    """
    Applies spin raising operator S+ to sites j and k, i.e. acts S+(j) S+(k) on a state.

    Parameters
    ----------
    state : String
        A string representing a state in ternary.
    j : Int
        An integer in the range [0, L-1] (where L is the total number of lattice sites) denoting a site j to raise.
    k : Int
        An integer in the range [0, L-1] (where L is the total number of lattice sites) denoting a site k to raise.

    Returns
    -------
    new_state : String
        The ternary representation of the new state after the spin raising operator S+ has acted on sites j and k.
    """
    
    state_rev = reverse(state)
    list_sites = collect(state_rev)

    list_sites[j + 1] = string(parse(Int, list_sites[j + 1]) - 1)
    list_sites[k + 1] = string(parse(Int, list_sites[k + 1]) - 1)

    new_state = join(reverse(list_sites))

    return new_state
end

function SzSz(state::String, j::Int, k::Int)
    """
    Applies Sz operator to sites j and k, and returns the coefficient multiplying the resultant (identical) state.

    Parameters
    ----------
    state : String
        A string representing a state in ternary.
    j : Int
        An integer in the range [0, L-1] (where L is the total number of lattice sites) denoting a site j where we apply the Sz operator.
    k : Int
        An integer in the range [0, L-1] (where L is the total number of lattice sites) denoting a site k where we apply the Sz operator.

    Returns
    -------
    coeff : Float64
        The coefficient resulting from the Sz operator acting on sites j and k.
    """
    
    state_rev = reverse(state)

    sj = parse(Int, state_rev[j + 1])
    sk = parse(Int, state_rev[k + 1])

    ms_values = [1, 0, -1]  # Spin z-components (possible m_s values) for spin-1

    coeff_j = ms_values[sj + 1]  # Adjusting for 1-based indexing
    coeff_k = ms_values[sk + 1]  # Adjusting for 1-based indexing

    coeff = coeff_j * coeff_k

    return coeff
end

function Ss(state::String, i::Int)
    """
    Applies squared spin operator (Sx + Sy + Sz)^2 to a site j, and returns the resultant states and their coefficients.

    Parameters
    ----------
    state : String
        A string representing a state in ternary.
    i : Int
        An integer in the range [0, L-1] (where L is the total number of lattice sites) denoting a site i where we apply the operator (Sx + Sy + Sz)^2.

    Returns
    -------
    result::Vector{Tuple{String, Complex{Float64}}}
        A vector of tuples of the form (state, coefficient), which result from applying the operator (Sx + Sy + Sz)^2 to our initial state at site i.
    """
    
    rev_state = reverse(state)     # Reverse string so lattice site 0 at beginning
    list_revstate = collect(rev_state)

    spin_val = rev_state[i + 1]  # Adjusting for 1-based indexing
    insq2 = 1 / sqrt(2)

    bprime1 = copy(list_revstate)
    bprime2 = copy(list_revstate)
    bprime3 = copy(list_revstate)

    bprime1_options = ["1", "0", "0"]
    bprime2_options = ["2", "2", "1"]

    coeff1_options = [insq2 * (1 + 1.0im), insq2 * (1 - 1.0im), -1.0im]
    coeff2_options = [1.0im, insq2 * (-1 - 1.0im), insq2 * (-1 + 1.0im)]

    bprime1[i + 1] = bprime1_options[parse(Int, spin_val) + 1]  # Adjusting for 1-based indexing
    bprime2[i + 1] = bprime2_options[parse(Int, spin_val) + 1]  # Adjusting for 1-based indexing

    coeff1 = coeff1_options[parse(Int, spin_val) + 1]  # Adjusting for 1-based indexing
    coeff2 = coeff2_options[parse(Int, spin_val) + 1]  # Adjusting for 1-based indexing

    # Rejoin lists to strings, and reverse to get in proper ternary form
    bprime1 = join(reverse(bprime1))
    bprime2 = join(reverse(bprime2))
    bprime3 = join(reverse(bprime3))

    coeff3 = 2

    return [(bprime1, coeff1), (bprime2, coeff2), (bprime3, coeff3)]
end

function construct_Hxx(L::Int, x_neighbors::Vector{Tuple{Int, Int}}, kept_ints::Vector{Int}, state_map::Vector{Int}, n_unique_list::Vector{Int}, filename::String, as_csr::Bool=true)
    """
    Constructs Hxx, i.e. the component of the Kitaev Hamiltonian which sums over all x-direction bonds of the
    honeycomb lattice.

    Parameters
    ----------
    L : Int
        The total number of lattice sites.
    x_neighbors : Vector{Tuple{Int, Int}}
        A list of tuples of site pairs which each constitute an x-bond.
    kept_ints : Vector{Int}
        An array of the base-10 integers corresponding to the representative states.
    state_map : Vector{Int}
        An array containing the representative state (in base-10) for each of the 3^L possible states.
    n_unique_list : Vector{Int}
        An array containing the number of unique mirror states for each representative state.
    filename : String
        Filename for the Hxx matrix which will be saved as a .npz file.
    as_csr : Bool, optional
        If True (default), the Hxx sparse matrix is converted from List of Lists (LIL) format to Compressed Sparse Row
        (CSR) format before saving.

    Returns
    -------
    Hxx : SparseMatrixCSC{Float64}
        The matrix Hxx in Compressed Sparse Row (CSR) format.
    """
    
    ndim = length(kept_ints)
    Hxx = sparsezeros(ndim, ndim)

    println("Constructing Hxx...\n")

    for i in 1:ndim
        a = ternary_pad(kept_ints[i], L)

        for (j, k) in x_neighbors
            bprime1 = SpSp(a, j, k)
            if !contains(bprime1, '3') && !contains(bprime1, '-')
                row = state_map[tern_to_base10(bprime1) + 1]  # Adjust for 1-based indexing
                col = i

                Ra = n_unique_list[i]
                Rb = n_unique_list[row]

                Hxx[row + 1, col + 1] += 0.25 * 2 * sqrt(Ra / Rb)
            end

            bprime2 = SpSm(a, j, k)
            if !contains(bprime2, '3') && !contains(bprime2, '-')
                row = state_map[tern_to_base10(bprime2) + 1]
                col = i

                Ra = n_unique_list[i]
                Rb = n_unique_list[row]

                Hxx[row + 1, col + 1] += 0.25 * 2 * sqrt(Ra / Rb)
            end

            bprime3 = SmSp(a, j, k)
            if !contains(bprime3, '3') && !contains(bprime3, '-')
                row = state_map[tern_to_base10(bprime3) + 1]
                col = i

                Ra = n_unique_list[i]
                Rb = n_unique_list[row]

                Hxx[row + 1, col + 1] += 0.25 * 2 * sqrt(Ra / Rb)
            end

            bprime4 = SmSm(a, j, k)
            if !contains(bprime4, '3') && !contains(bprime4, '-')
                row = state_map[tern_to_base10(bprime4) + 1]
                col = i

                Ra = n_unique_list[i]
                Rb = n_unique_list[row]

                Hxx[row + 1, col + 1] += 0.25 * 2 * sqrt(Ra / Rb)
            end
        end
    end

    if as_csr
        println("Converting Hxx to csr format...\n")
        Hxx = sparse(Hxx)
    end

    println("Saving Hxx...\n")
    save_npz(filename, Hxx)

    return Hxx
end

function construct_Hyy(L::Int, y_neighbors::Vector{Tuple{Int, Int}}, kept_ints::Vector{Int}, state_map::Vector{Int}, n_unique_list::Vector{Int}, filename::String, as_csr::Bool=true)
    """
    Constructs Hyy, i.e. the component of the Kitaev Hamiltonian which sums over all y-direction bonds of the
    honeycomb lattice.
    
    Parameters
    ----------
    L : Int
        The total number of lattice sites.
    y_neighbors : Vector{Tuple{Int, Int}}
        A list of tuples of site pairs which each constitute a y-bond.
    kept_ints : Vector{Int}
        An array of the base-10 integers corresponding to the representative states.
    state_map : Vector{Int}
        An array containing the representative state (in base-10) for each of the 3^L possible states.
    n_unique_list : Vector{Int}
        An array containing the number of unique mirror states for each representative state.
    filename : String
        Filename for the Hyy matrix which will be saved as a .npz file.
    as_csr : Bool, optional
        If True (default), the Hyy sparse matrix is converted from List of Lists (LIL) format to Compressed Sparse Row
        (CSR) format before saving.

    Returns
    -------
    Hyy : SparseMatrixCSC{Float64}
        The matrix Hyy in Compressed Sparse Row (CSR) format.
    """
    
    ndim = length(kept_ints)
    Hyy = sparsezeros(ndim, ndim)

    println("Constructing Hyy...\n")

    for i in 1:ndim
        a = ternary_pad(kept_ints[i], L)

        for (j, k) in y_neighbors
            bprime1 = SpSp(a, j, k)
            if !contains(bprime1, '3') && !contains(bprime1, '-')
                row = state_map[tern_to_base10(bprime1) + 1]
                col = i

                Ra = n_unique_list[i]
                Rb = n_unique_list[row]

                Hyy[row + 1, col + 1] -= 0.25 * 2 * sqrt(Ra / Rb)
            end

            bprime2 = SpSm(a, j, k)
            if !contains(bprime2, '3') && !contains(bprime2, '-')
                row = state_map[tern_to_base10(bprime2) + 1]
                col = i

                Ra = n_unique_list[i]
                Rb = n_unique_list[row]

                Hyy[row + 1, col + 1] += 0.25 * 2 * sqrt(Ra / Rb)
            end

            bprime3 = SmSp(a, j, k)
            if !contains(bprime3, '3') && !contains(bprime3, '-')
                row = state_map[tern_to_base10(bprime3) + 1]
                col = i

                Ra = n_unique_list[i]
                Rb = n_unique_list[row]

                Hyy[row + 1, col + 1] += 0.25 * 2 * sqrt(Ra / Rb)
            end

            bprime4 = SmSm(a, j, k)
            if !contains(bprime4, '3') && !contains(bprime4, '-')
                row = state_map[tern_to_base10(bprime4) + 1]
                col = i

                Ra = n_unique_list[i]
                Rb = n_unique_list[row]

                Hyy[row + 1, col + 1] -= 0.25 * 2 * sqrt(Ra / Rb)
            end
        end
    end

    if as_csr
        println("Converting Hyy to csr format...\n")
        Hyy = sparse(Hyy)
    end

    println("Saving Hyy...\n")
    save_npz(filename, Hyy)

    return Hyy
end

function construct_Hzz(L::Int, z_neighbors::Vector{Tuple{Int, Int}}, kept_ints::Vector{Int}, filename::String, as_csr::Bool=true)
    """
    Constructs Hzz, i.e. the component of the Kitaev Hamiltonian which sums over all z-direction bonds of the
    honeycomb lattice.

    Parameters
    ----------
    L : Int
        The total number of lattice sites.
    z_neighbors : Vector{Tuple{Int, Int}}
        A list of tuples of site pairs which each constitute a z-bond.
    kept_ints : Vector{Int}
        An array of the base-10 integers corresponding to the representative states.
    filename : String
        Filename for the Hzz matrix which will be saved as a .npz file.
    as_csr : Bool, optional
        If True (default), the Hzz sparse matrix is converted from List of Lists (LIL) format to Compressed Sparse Row
        (CSR) format before saving.

    Returns
    -------
    Hzz : SparseMatrixCSC{Float64}
        The matrix Hzz in Compressed Sparse Row (CSR) format.
    """

    ndim = length(kept_ints)
    Hzz = sparsezeros(ndim, ndim)
    println("Constructing Hzz...\n")
    for i in 1:ndim
        a = ternary_pad(kept_ints[i], L)
        for (j, k) in z_neighbors
            coeff = SzSz(a, j, k)
            Hzz[i, i] += coeff
        end
    end
    if as_csr
        println("Converting Hzz to csr format...\n")
        Hzz = sparse(Hzz)
    end
    println("Saving Hzz...\n")
    save_npz(filename, Hzz)
    return Hzz
end

function construct_HD(L::Int, kept_ints::Vector{Int}, state_map::Vector{Int}, n_unique_list::Vector{Int}, filename::String, as_csr::Bool=true)
    """
    Constructs HD, i.e. the component of the Hamiltonian which describes a single-ion anisotropy in the [1, 1, 1]
    direction. The HD operator involves a sum over all sites of D * (Sx + Sy + Sz)^2, where D is the magnitude of the
    single-ion anisotropy.

    Parameters
    ----------
    L : Int
        The total number of lattice sites.
    kept_ints : Vector{Int}
        An array of the base-10 integers corresponding to the representative states.
    state_map : Vector{Int}
        An array containing the representative state (in base-10) for each of the 3^L possible states.
    n_unique_list : Vector{Int}
        An array containing the number of unique mirror states for each representative state.
    filename : String
        Filename for the HD matrix which will be saved as a .npz file.
    as_csr : Bool, optional
        If True (default), the HD sparse matrix is converted from List of Lists (LIL) format to Compressed Sparse Row
        (CSR) format before saving.

    Returns
    -------
    HD : SparseMatrixCSC{Complex{Float64}}
        The matrix HD in Compressed Sparse Row (CSR) format.
    """

    ndim = length(kept_ints)
    HD = sparsezeros(ndim, ndim, Complex{Float64})

    println("Constructing HD...\n")

    for i in 1:ndim
        a = ternary_pad(kept_ints[i], L)

        for l in 1:L
            Ss_result = Ss(a, l)

            bprime1, coeff1 = Ss_result[1]
            col = i
            row = state_map[tern_to_base10(bprime1) + 1]

            Ra = n_unique_list[i]
            Rb = n_unique_list[row]

            HD[row + 1, col + 1] += coeff1 * sqrt(Ra / Rb)

            bprime2, coeff2 = Ss_result[2]
            row = state_map[tern_to_base10(bprime2) + 1]

            Rb = n_unique_list[row];

            HD[row + 1, col + 1] += coeff2 * sqrt(Ra / Rb)

            # bprime3 = a, so row = col, and coeff3 = 2
            HD[col + 1, col + 1] += 2
        end
    end

    if as_csr
        println("Converting HD to csr format...\n")
        HD = sparse(HD)
    end

    println("Saving HD...\n")
    save_npz(filename, HD)

    return HD
end

function construct_hamiltonian(Kx::Float64, Ky::Float64, Kz::Float64, D::Float64, Hxx_file::String, Hyy_file::String, Hzz_file::String, HD_file::String)
    """
    Loads from file separate components of the Hamiltonian (Hxx, Hyy, Hzz, HD) and constructs the full Hamiltonian
    matrix.

    Parameters
    ----------
    Kx : Float64
        The Kitaev coupling along x-bonds.
    Ky : Float64
        The Kitaev coupling along y-bonds.
    Kz : Float64
        The Kitaev coupling along z-bonds.
    D : Float64
        The magnitude of the single-ion anisotropy term, i.e. the coefficient which multiplies HD.
    Hxx_file: String
        Saved file (in .npz format) for the Hxx matrix.
    Hyy_file: String
        Saved file (in .npz format) for the Hyy matrix.
    Hzz_file: String
        Saved file (in .npz format) for the Hzz matrix.
    HD_file: String
        Saved file (in .npz format) for the HD matrix.

    Returns
    -------
    H : SparseMatrixCSC{Complex{Float64}}
        The full sparse matrix Kitaev Hamiltonian corresponding to H = (Kx * Hxx) + (Ky * Hyy) + (Kz * Hzz) + (D * HD).
    """

    println("Loading Hxx...\n")
    Hxx = load_npz(Hxx_file)

    println("Loading Hyy...\n")
    Hyy = load_npz(Hyy_file)

    println("Loading Hzz...\n")
    Hzz = load_npz(Hzz_file)

    println("Loading HD...\n")
    HD = load_npz(HD_file)

    println("Constructing full Hamiltonian...\n")
    H = (Kx * Hxx) + (Ky * Hyy) + (Kz * Hzz) + (D * HD)

    println("Hamiltonian constructed\n")

    return H
end
```

### Key Changes and Considerations:

1. **Function Definitions**: Each function is defined with clear parameters and return types, using Julia's type annotations.
  
2. **1-Based Indexing**: Julia uses 1-based indexing, so adjustments have been made where necessary (e.g., `row + 1`, `col + 1`).

3. **Sparse Matrices**: The functions utilize Julia's `SparseArrays` for creating and manipulating sparse matrices. The `sparsezeros` function initializes a sparse matrix with zeros.

4. **String Operations**: The use of `contains()` checks if a string contains certain characters, which is a direct translation from Python.

5. **Complex Numbers**: The matrices that require complex numbers use `Complex{Float64}` to denote their type.

6. **Saving and Loading**: The functions `save_npz` and `load_npz` are assumed to be defined for handling `.npz` file operations. You'll need to implement or import these functions based on how you want to manage file operations in Julia.

### Usage Example

You would use the defined functions as follows:

```
