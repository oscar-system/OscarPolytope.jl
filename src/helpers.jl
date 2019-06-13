matrix_for_polymake(x::Nemo.fmpz_mat) = Matrix{BigInt}(x)
matrix_for_polymake(x::Nemo.fmpq_mat) = Matrix{Rational{BigInt}}(x)
matrix_for_polymake(x::AbstractMatrix{<:Integer}) = x
matrix_for_polymake(x::AbstractMatrix{<:Rational{<:Integer}}) = x

augment(vec::AbstractVector{T}, val) where T = vcat(T(val), vec)
augment(mat::AbstractMatrix, vec::AbstractVector) = hcat(vec, mat)

homogenize(::Type{T}, vec::AbstractVector, val::T=T(0)) where T = augment(T.(vec), val)
homogenize(::Type{T}, mat::AbstractMatrix, val::T=T(1)) where T = augment(T.(mat), repeat([val], size(mat,1)))
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
