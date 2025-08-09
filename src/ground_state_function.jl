# Import necessary functions
# Assuming equivalent functions are defined in the appropriate module or defined here

# Define the functions required for the calculations

function ground_state(H::AbstractMatrix, L::Int, k::Int=1, ncv::Int=20)
    """
    Finds the ground state psi_0 and the ground state energy E_0, by applying
    the Lanczos algorithm to the matrix Hamiltonian H.

    Parameters
    ----------
    H : AbstractMatrix
        The full sparse matrix Hamiltonian.
    L : Int
        The total number of lattice sites. Used to calculate the energy per site.
    k : Int, optional
        Number of eigenvalues and eigenvectors desired. By default k=1 to find state only.
    ncv : Int, optional
        Number of Lanczos vectors generated.

    Returns
    -------
    E_0 : Float64
        The ground state energy per site.
    psi_0 : Vector{Complex{Float64}}
        The ground state, i.e. a 1D array with complex coefficients.
    """
    
    E, V = eigen(H; k=k, which=:smallestabs)

    # Get ground state energy per site (i.e. divide by L) and corresponding eigenvector
    E_0 = E[1] / L
    psi_0 = V[:, 1]

    println("Ground state energy (per site): ", E_0)

    return E_0, psi_0
end

function s_squared(HD::AbstractMatrix, L::Int, psi_0::Vector{Complex{Float64}})
    """
    Finds the expectation value of (Sx + Sy + Sz)^2 in the ground state (per site).

    Parameters
    ----------
    HD : AbstractMatrix
        The sparse matrix HD corresponding to the single-ion anisotropy term in the Hamiltonian.
    L : Int
        The total number of lattice sites. Used to calculate <(Sx + Sy + Sz)^2> per site.
    psi_0 : Vector{Complex{Float64}}
        The ground state, i.e. a 1D array with complex coefficients.

    Returns
    -------
    Ss: Float64
        The expectation value of (Sx + Sy + Sz)^2 in the ground state (per site).
    """
    
    Ss = conj(psi_0)' * HD * psi_0

    Ss = real(Ss) / L  # Value per site

    return Ss
end

function entanglement_entropy(L::Int, psi::Vector{Complex{Float64}}, state_map::Vector{Int}, n_unique_list::Vector{Int})
    """
    Finds the bipartite entanglement entropy when the system is divided into two equal halves.
    Performs a singular value decomposition (SVD) to obtain the singular values.

    Parameters
    ----------
    L : Int
        The total number of lattice sites.
    psi : Vector{Complex{Float64}}
        The state of the system, i.e. a 1D array with complex coefficients. Typically the ground state psi_0.
    state_map : Vector{Int}
        An array containing the representative state (in base-10) for each of the 3^L possible states.
    n_unique_list : Vector{Int}
        An array containing the number of unique mirror states for each representative state.

    Returns
    -------
    entropy: Float64
        The bipartite entanglement entropy when the system is divided into two equal halves.
    """
    
    M_dim = 3^(L ÷ 2)
    M = zeros(Complex{Float64}, M_dim, M_dim)

    println("Constructing $M_dim x $M_dim matrix for SVD...\n")

    for i in 1:M_dim
        for j in 1:M_dim
            state_left = ternary_pad(i - 1, L ÷ 2)  # Adjusting for 1-based indexing
            state_right = ternary_pad(j - 1, L ÷ 2) # Adjusting for 1-based indexing

            state = state_left .+ state_right
            state_int = tern_to_base10(state)

            rep_state_index = state_map[state_int + 1]  # Adjusting for 1-based indexing
            ci = psi[rep_state_index + 1]  # Adjusting for 1-based indexing
            ri = n_unique_list[rep_state_index + 1]  # Adjusting for 1-based indexing

            M[i, j] = ci / sqrt(ri)
        end
    end

    println("Performing SVD...\n")

    # Get array of singular values:
    S = svd(M).S

    println("Calculating entanglement entropy...\n")

    entropy = 0.0
    for value in S
        entropy += -1 * (value^2) * log(value^2)
    end

    println("Entropy calculation complete \n")

    return entropy
end
```

### Explanation of Key Changes and Considerations:

1. **Function Signatures**: In Julia, we define function signatures with type annotations, specifying the types of inputs and outputs.

2. **Matrix Operations**: The operations involving matrices are adjusted to use Julia's matrix multiplication syntax (`*` for matrix multiplication, `.*` for element-wise multiplication).

3. **Complex Numbers**: Julia has built-in support for complex numbers, so we use `Complex{Float64}` for complex coefficients.

4. **SVD**: The singular value decomposition is performed using `svd()` in Julia, which returns a structure containing singular values.

5. **1-based Indexing**: Julia arrays are 1-based, so all indices are adjusted accordingly (for example, `i - 1` when accessing an array).

6. **Print Statements**: The `println` function is used for output, similar to Python's `print`.

7. **Assumed Functions**: The functions `ternary_pad` and `tern_to_base10` are assumed to exist and perform the same operations as in the original Python code.

### Usage Example

You would use the defined functions as follows:

```