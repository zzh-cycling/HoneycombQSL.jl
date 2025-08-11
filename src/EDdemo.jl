# Hamiltonian is defined as H = -∑_i,j (X_i X_j + Y_i Y_j + Z_i Z_j), on a honeycomb lattice with periodic boundary conditions, where X, Y, Z are Pauli matrices acting on the spins at sites i and j. Here we use the parallelogram (or skewed-rectangle) cluster, tilted square cluster, rhombic cluster lattice arranging way, whatever you call it, to define the Hamiltonian.  Such arranging is easy to define the K-space, utilizing the translation symmetry of the honeycomb lattice.

# Another way is using hexagonal (or honeycomb) cluster, diamond-shaped cluster, zig-zag cluster, which is born to be OBC. If PBC, we need to do twisted/skewed PBC or modular boundary condition.

function honeycomb_strings(m::Int, n::Int, pbc::Bool=true; cluster::Symbol= :rhombic)
	@assert cluster in [:rhombic, :tilted_square, :parallelogram, :hexagonal, :zigzag]
	li = reshape(LinearIndices((m, n)),(n,m))'
	xstrings = Vector{Int}[]
	ystrings = Vector{Int}[]
	zstrings = Vector{Int}[]

	if cluster ∈ [:rhombic, :tilted_square, :parallelogram]
		# Rhombic cluster, tilted square cluster, parallelogram cluster
		for i=1:m, j=1:n
			label = li[i, j]
			push!(ystrings, [2*label-2,2*label-1])
		end
		for i=1:m-1, j=1:n
			label1= li[i, j]
			label2= li[i+1, j]
			push!(xstrings, [2*label1-1, 2*label2-2])	
		end
		for i=1:m, j=1:n-1
			label1 = li[i, j]
			label2 = li[i, j+1]
			push!(zstrings, [2*label1-1, 2*label2-2])
		end
	
		if pbc
			for j=1:n
				label1= li[m, j]
				label2 = li[1, j]
				push!(xstrings, [2*label1-1, 2*label2-2])
			end
	
			for i=1:m
				label1= li[i, 1]
				label2 = li[i, n]
				push!(zstrings, [2*label1-2, 2*label2-1])
			end
		end
	else

	end

	return xstrings, ystrings, zstrings
end

function kagome_strings(m::Int, n::Int, pbc::Bool=true; cluster::Symbol= :rhombic)
	# Kagome lattice strings, not implemented yet.
	error("Kagome lattice strings not implemented yet.")
end

⊗(A::AbstractArray, B::AbstractArray) = kron(A, B)

function kitaev_hamiltonian(m::Int, n::Int)
	xstrings, ystrings, zstrings = honeycomb_strings(m, n)
	X= [0 1; 1 0]
	Y= [0 -1; 1 0] # Note not need to substract YY
	Z= [1 0; 0 -1]
	Id= [1 0; 0 1]
	
	l=2*m*n
	H=zeros(Int64, 2^l, 2^l)
	for strings in xstrings
		cup=fill(Id, l)
		cup[(strings .+1)] = fill(X, 2)
		H-=foldr(⊗,cup)
	end

	for strings in ystrings
		cup=fill(Id, l)
		cup[(strings .+1)] = fill(Y, 2)
		H+=foldr(⊗,cup)
	end

	for strings in zstrings
		cup=fill(Id, l)
		cup[(strings .+1)] = fill(Z, 2)
		H-=foldr(⊗,cup)
	end
	return H
end

function kitaev_hamiltonian_sparse(m::Int, n::Int)
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
		H-=foldr(⊗,cup)
    end

    for strings in ystrings
        cup=fill(Id, l)
		cup[(strings .+1)] = fill(Y, 2)
		H+=foldr(⊗,cup)
    end

    for strings in zstrings
        cup=fill(Id, l)
		cup[(strings .+1)] = fill(Z, 2)
		H-=foldr(⊗,cup)
    end
	
    return H
end

function Heisenberg_hamiltonian_sparse(m::Int, n::Int; lattice::Symbol=:honeycomb)
	# Heisenberg Hamiltonian for kagome and honeycomb lattice
	if lattice == :honeycomb
		xstrings, ystrings, zstrings = honeycomb_strings(m, n)
	elseif lattice == :kagome
		xstrings, ystrings, zstrings = kagome_strings(m, n)
	end

	X = sparse([0 1; 1 0])
	Y = sparse([0 -1; 1 0])
	Z = sparse([1 0; 0 -1])
	Id = sparse([1 0; 0 1])

	l = 2*m*n
	H = spzeros(2^l, 2^l)  

	for strings in xstrings
		cup=fill(Id, l)
		cup[(strings .+1)] = fill(X, 2)
		H+=foldr(⊗,cup)
	end

	for strings in ystrings
		cup=fill(Id, l)
		cup[(strings .+1)] = fill(Y, 2)
		H+=foldr(⊗,cup)
	end

	for strings in zstrings
		cup=fill(Id, l)
		cup[(strings .+1)] = fill(Z, 2)
		H+=foldr(⊗,cup)
	end

	return H
end

function loop_path(i::Int, j::Int, m::Int, n::Int, pbc::Bool=true)
	# First for PBC. m is the number of rows, n is the number of columns, i,j is the plaquette position (or cell.)
	@assert(i in 1:m && j in 1:n)
	li = reshape(LinearIndices((m, n)),(n,m))'

	modi = mod1(i+1, m)
	modj = mod1(j+1, n)
	if i<m || j < n 
		path = [2*li[i,j]-1, 2*li[i,modj]-2, 2*li[i,modj]-1, 
		2*li[modi,j]-2, 2*li[modi,j]-1, 2*li[modi,j]] 
	end
	if pbc && (i==m || j==n)
		if i==m && j<n
			path = [2*li[i,j]-1, 2*li[i,modj]-2, 2*li[i,modj]-1, 
			2*li[1,j]-2, 2*li[1,j]-1, 2*li[1,j]] 
		elseif i<m && j==n
			path = [2*li[i,j]-1, 2*li[i,modj]-2, 2*li[i,modj]-1, 
			2*li[i+1,j]-2, 2*li[i+1,j]-1, 2*li[i+1,1]-2] 
		else
			path = [2*li[m,j]-1, 2*li[m,modj]-2, 2*li[m,modj]-1,
			2*li[1,j]-2, 2*li[1,j]-1, 2*li[1,1]-2] 
		end
	end

	return path
end

function wilson12(m::Int, n::Int)
	li = reshape(LinearIndices((m, n)),(n,m))'
	l=2*m*n
	X = sparse([0 1; 1 0])
	Y = sparse([0 -1; 1 0])  
	Z = sparse([1 0; 0 -1]) 
	Id = sparse([1 0; 0 1])
	cup1=fill(Id, l)
	cup2=fill(Id, l)

	path1 = foldl(vcat, [[2*li[i,1]-2, 2*li[i, 1]-1] for i in 1:m]) .+1
	path2 = foldl(vcat, [[2*li[1,j]-2, 2*li[1,j]-1] for j in 1:n]) .+1
	cup1[path1] = vcat([X], (-1)^(m-1)*fill(Z, 2*m-2), [X])
	cup2[path2] = vcat([X], (-1)^(n-1)*fill(Y, 2*n-2), [X])

	wilson1=foldr(⊗,cup1)
	wilson2=foldr(⊗,cup2)
	return wilson1, wilson2
end

function loop_op(i::Int, j::Int, m::Int, n::Int)
	# First for PBC. m is the number of rows, n is the number of columns, i,j is the plaquette position (or cell.).
	@assert(i in 1:m && j in 1:n)
	
	l=2*m*n
	X = sparse([0 1; 1 0])
    Y = sparse([0 -1; 1 0])  
    Z = sparse([1 0; 0 -1]) 
    Id = sparse([1 0; 0 1])
	cup=fill(Id, l)

	# Here we use YY definition, which is different from Y, has a additional minus sign. W = -Y1 X2 Z3 Y10 X9 Z8. In kitaev_113 order, X1,Y2,Z3,Y4,X5,Z6.
	path = loop_path(i, j, m, n)
	cup[(path .+1)] = [-Y, X, Z, Z, X, Y]

	loop = foldr(⊗,cup)
	return loop
end

