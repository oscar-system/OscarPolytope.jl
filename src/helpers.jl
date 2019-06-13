matrix_for_polymake(x::Nemo.fmpz_mat) = Matrix{BigInt}(x)
matrix_for_polymake(x::Nemo.fmpq_mat) = Matrix{Rational{BigInt}}(x)
matrix_for_polymake(x::AbstractMatrix{<:Integer}) = x
matrix_for_polymake(x::AbstractMatrix{<:Rational{<:Integer}}) = x
matrix_for_polymake(x::Polymake.pm_MatrixAllocated{Polymake.pm_Rational}) = x

augment(vec::V, val) where {T, V<:AbstractVector{T}} = V(vcat(val, vec))
augment(mat::M, vec::AbstractVector) where {T, M<:AbstractMatrix{T}} = M(hcat(vec, mat))

homogenize(::Type{T}, vec::AbstractVector, val=T(0)) where T = augment(T.(vec), val)
homogenize(::Type{T}, mat::AbstractMatrix, val=T(1)) where T = augment(T.(mat), fill(val, size(mat,1)))
homogenize(mat::AbstractVecOrMat{T}, val=T(1)) where T = homogenize(T, mat,
val)

dehomogenize(vec::AbstractVector) = vec[2:end]
dehomogenize(mat::AbstractMatrix) = mat[:, 2:end]

function property_is_computed(P::Polymake.pm_perl_ObjectAllocated, S::Symbol)
   pv = Polymake.internal_call_method("lookup", P, Any[string(S)])
   return nothing != Polymake.convert_from_property_value(pv)
end
function property_is_computed(HP::HomogeneousPolyhedron, S::Symbol)
   return property_is_computed(HP.polymakePolytope, S)
end
function property_is_computed(P::Polyhedron, S::Symbol)
   return property_is_computed(P.homogeneous_polyhedron, S)
end

function decompose_vdata(A)
   ncols = (size(A))[2]
   vertexIndices = Int[]
   rayIndices = Int[]
   for i=1:ncols
      if A[1, i] == 0
         push!(rayIndices, i)
      elseif A[1, i] == 1
         push!(vertexIndices, i)
      end
   end
   return (A[2:end, vertexIndices], A[2:end, rayIndices])
end

function decompose_hdata(A)
   (-A[:,2:end], A[:,1:1])
end
