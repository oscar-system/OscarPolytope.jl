function augment(vec::AbstractVector, val)
   s = size(vec)
   res = similar(vec, (s[1]+1,))
   res[1] = val
   res[2:end] = vec
   return res
end

function augment(mat::AbstractMatrix, vec::AbstractVector)
   s = size(mat)
   res = similar(mat, (s[1], s[2]+1))
   res[:, 1] = vec
   res[:, 2:end] = mat
   return res
end

homogenize(vec::AbstractVector, val::Number=0) = augment(vec, val)
homogenize(mat::AbstractMatrix, val::Number=1) = augment(mat, fill(val, size(mat,1)))

dehomogenize(vec::AbstractVector) = vec[2:end]
dehomogenize(mat::AbstractMatrix, rowidcs=1:size(mat,1)) = mat[rowidcs, 2:end]

iscomputed(HP::HomogeneousPolyhedron, property::Symbol) = Polymake.exists(HP.polymakePolytope, string(property))
iscomputed(P::Polyhedron, property::Symbol) = iscomputed(P.homogeneous_polyhedron, property)

#assuming polymake, row-major storage
vertex_indices(A::AbstractMatrix) = filter(i->isone(A[i,1]), 1:size(A,1))
ray_indices(A::AbstractMatrix) = filter(i->iszero(A[i,1]), 1:size(A,1))
vertices(A::AbstractMatrix) = dehomogenize(A, vertex_indices(A))
rays(A::AbstractMatrix) = dehomogenize(A, ray_indices(A))
decompose_hdata(A) = (-dehomogenize(A), A[:,1])

function Base.convert(::Type{Polymake.pm_Integer},
   x::Union{Nemo.fmpz, Nemo.fmpq})
   return Polymake.pm_Integer(BigInt(x))
end

function Base.convert(::Type{Polymake.pm_Rational},
   x::Union{Nemo.fmpz, Nemo.fmpq})
   return Polymake.pm_Rational(convert(Rational{BigInt}, x))
end

function Base.convert(::Type{Polymake.PolymakeType},
   x::Nemo.MatrixElem{<:Union{Nemo.RingElem, Integer}})
   return Polymake.pm_Matrix{Polymake.pm_Integer}(Matrix(x))
end

function Base.convert(::Type{Polymake.PolymakeType},
   x::Nemo.MatrixElem{<:Union{Nemo.FracElem{Nemo.fmpz}, Rational}})
   return Polymake.pm_Matrix{pm_Rational}(x)
end
