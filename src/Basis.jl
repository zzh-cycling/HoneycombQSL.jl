function honeycomb_basis(m::Int,n::Int, Wplis::Vector{Int}=Int[]; spin::Float64=0.5)
    # Generate the basis for the honeycomb lattice Kitaev model for spin spin.
    # m: number of rows, n: number of columns
    l = 2 * m * n  # total number of spins
    d = round(Int, 2*spin+1)  # dof
    if isempty(Wplis)
        return [DitStr{d, l}(i) for i in 0:2^l-1]
    else
        return [DitStr{d, l}(i) for i in 0:2^l-1 if all(x -> x in Wplis, BitBasis.bitstring(i, l))]
    end
end

function kagome_basis(m::Int, n::Int, Wplis::Vector{Int}=Int[]; spin::Float64=0.5)
    # Generate the basis for the kagome lattice Kitaev model for spin spin.
    # m: number of rows, n: number of columns
    l = 2 * m * n  # total number of spins
    d = round(Int, 2*spin+1)  # dof
    if isempty(Wplis)
        return [DitStr{d, l}(i) for i in 0:2^l-1]
    else
        return [DitStr{d, l}(i) for i in 0:2^l-1 if all(x -> x in Wplis, BitBasis.bitstring(i, l))]
    end
    
end

function honeycomb_basis(m::Int,n::Int, Wplis::Vector{Int}, k::Vector{Int}, spin::Float64=0.5)
    # k is the momentum space basis.
    l = 2 * m * n  # total number of spins
    d = round(Int, 2*spin+1)  # dof
    return [DitStr{d, l}(i) for i in 0:2^l-1]
end

function kitaev_ham_map(m::Int, n::Int, string::Vector{Int}, stringtype::Symbol; spin::Float64=0.5)
    # Generate the mapping for the honeycomb lattice Kitaev model.
    # m: number of rows, n: number of columns, spin: spin value
    l = 2 * m * n  # total number of spins
    d = round(Int, 2*spin+1)  # dof
    return Dict{Int, Int}(i => i for i in 0:2^l-1)
end

function loop_map(i::Int, j::Int, m::Int, n::Int; spin::Float64=0.5)
    # Generate the mapping for the loop operators.
    # m: number of rows, n: number of columns, i,j: plaquette position
    l = 2 * m * n  # total number of spins
    d = round(Int, 2*spin+1)  # dof
    basis = honeycomb_basis(m, n; spin=spin)
    path = loop_path(i, j, m, n)
    new_state = [copy(b) for b in basis]
    return new_state
end

function kitaev_honeycomb_ham(m::Int, n::Int, Jx::Float64=1.0, Jy::Float64=1.0, Jz::Float64=1.0; h::Float64=0.0, spin::Float64=0.5, cluster=:triangular)
    # Generate the Hamiltonian for the honeycomb lattice Kitaev model.
    # m: number of rows, n: number of columns, spin: spin value, Jx, Jy, Jz: coupling constants, h: external magnetic field, cluster: type of cluster
    xstrings, ystrings, zstrings = honeycomb_strings(m, n)
    
    X = sparse([0 1; 1 0])
    Y = sparse([0 -1; 1 0])  
    Z = sparse([1 0; 0 -1]) 
    Id = sparse([1 0; 0 1])
    
    l = 2*m*n
    H = spzeros(2^l, 2^l)  # 使用复数矩阵

    for strings in xstrings
        cup=fill(Id, l)
        cup[(strings .+1)] = fill(X, 2)
        H -= Jx * foldr(⊗,cup)
    end

    for strings in ystrings
        cup=fill(Id, l)
        cup[(strings .+1)] = fill(Y, 2)
        H += Jy * foldr(⊗,cup)
    end

    for strings in zstrings
        cup=fill(Id, l)
        cup[(strings .+1)] = fill(Z, 2)
        H -= Jz * foldr(⊗,cup)
    end
    
    return H
end