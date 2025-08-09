function honeycomb_strings(m::Int, n::Int, pbc::Bool=true)
	li = reshape(LinearIndices((m, n)),(n,m))'
	xstrings = Vector{Int}[]
	ystrings = Vector{Int}[]
	zstrings = Vector{Int}[]

	for i=1:m, j=1:n
		label = li[i, j]
		push!(xstrings, [2*label-1,2*label])
	end
	for i=1:m-1, j=1:n
		label1= li[i, j]
		label2= li[i+1, j]
		push!(ystrings, [2*label1, 2*label2-1])	
	end
	for i=1:m, j=1:n-1
		label1 = li[i, j]
		label2 = li[i, j+1]
		push!(zstrings, [2*label1, 2*label2-1])
	end

	if pbc
		for j=1:n
			label1= li[m, j]
			label2 = li[1, j]
			push!(ystrings, [2*label1, 2*label2-1])
		end

		for i=1:m
			label1= li[i, 1]
			label2 = li[i, n]
			push!(zstrings, [2*label1-1, 2*label2])
		end
	end

	return xstrings, ystrings, zstrings
end

⊗(A::AbstractArray, B::AbstractArray) = kron(A, B)

function kitaev_hamiltonian(m::Int, n::Int)
	xstrings, ystrings, zstrings = honeycomb_strings(m, n)
	X= [0 1; 1 0]
	Y= [0 -1; 1 0] # Note need not substract YY
	Z= [1 0; 0 -1]
	Id= [1 0; 0 1]
	
	l=2*m*n
	H=zeros(Int64, 2^l, 2^l)
	for i in xstrings
		cup=fill(Id, l)
		for j in i
			cup[j]=X
		end
		H-=foldr(⊗,cup)
	end

	for i in ystrings
		cup=fill(Id, l)
		for j in i
			cup[j]=Y
		end
		H+=foldr(⊗,cup)
	end

	for i in zstrings
		cup=fill(Id, l)
		for j in i
			cup[j]=Z
		end
		H-=foldr(⊗,cup)
	end
	return H
end

function kitaev_hamiltonian_sparse(m::Int, n::Int)
    xstrings, ystrings, zstrings = honeycomb_strings(m, n)
    
    σx = sparse([0 1; 1 0])
    σy = sparse([0 -1; 1 0])  
    σz = sparse([1 0; 0 -1]) 
    Id = sparse([1 0; 0 1])
    
    l = 2*m*n
    H = spzeros(2^l, 2^l)  # 使用复数矩阵

    for sites in xstrings
        cup=fill(Id, l)
		for j in sites
			cup[j]=σx
		end
		H-=foldr(⊗,cup)
    end

    for sites in ystrings
        cup=fill(Id, l)
		for j in sites
			cup[j]=σy
		end
		H+=foldr(⊗,cup)
    end

    for sites in zstrings
        cup=fill(Id, l)
		for j in sites
			cup[j]=σz
		end
		H-=foldr(⊗,cup)
    end

    return H
end

function flux_path(i::Int, j::Int, m::Int, n::Int, pbc::Bool=true)
	@assert(i in 1:m && j in 1:n)
	li = reshape(LinearIndices((m, n)),(n,m))'
	l=2*m*n

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

	return path .+1
end

function Wilson12(i::Int, j::Int, m::Int, n::Int)
	# First for PBC. m is the number of rows, n is the number of columns, i,j is the plaquette position (or cell.)
	@assert(i in 1:m && j in 1:n)
	
	l=2*m*n
	σx = sparse([0 1; 1 0])
	σy = sparse([0 -1; 1 0])  
	σz = sparse([1 0; 0 -1]) 
	Id = sparse([1 0; 0 1])
	cup=fill(Id, l)

	path = flux_path(i, j, m, n)
	cup[path] = [σx, σy, σz, σz, σy, σx]

	wilson=foldr(⊗,cup)
	return wilson
end

function flux(i::Int, j::Int, m::Int, n::Int)
	# First for PBC. m is the number of rows, n is the number of columns, i,j is the plaquette position (or cell.)
	@assert(i in 1:m && j in 1:n)
	
	l=2*m*n
	σx = sparse([0 1; 1 0])
    σy = sparse([0 -1; 1 0])  
    σz = sparse([1 0; 0 -1]) 
    Id = sparse([1 0; 0 1])
	cup=fill(Id, l)

	path = flux_path(i, j, m, n)
	cup[path] = [σx, σy, σz, σz, σy, σx]

	flux=foldr(⊗,cup)
	return flux
end

function majorana(m::Int, n::Int)
	
end