# matrix_for_polymake(x::Nemo.fmpz_mat) = Matrix{BigInt}(x)
# matrix_for_polymake(x::Nemo.fmpq_mat) = Matrix{Rational{BigInt}}(x)
matrix_for_polymake(x::AbstractMatrix{<:Integer}) = x
matrix_for_polymake(x::AbstractMatrix{<:Rational{<:Integer}}) = x
matrix_for_polymake(x::Polymake.Matrix) = x
matrix_for_polymake(x::AbstractArray) = x

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
dehomogenize(mat::AbstractMatrix) = mat[:, 2:end]
#
function property_is_computed(P::Polymake.BigObjectAllocated, S::Symbol)
   pv = Polymake.internal_call_method("lookup", P, Any[string(S)])
   return nothing != Polymake.convert_from_property_value(pv)
end

"""
   decompose_vdata(A::AbstractMatrix)

Given a (homoegenuous) polymake matrix split into vertices and rays and dehomogenize.
"""
function decompose_vdata(A::AbstractMatrix)
   vertex_indices = findall(!iszero, view(A, :, 1))
   ray_indices = findall(iszero, view(A, :, 1))
   return (vertices = A[vertex_indices, 2:end], rays = A[ray_indices, 2:end])
end

function decompose_hdata(A)
   (A = -A[:,2:end], b = A[:,1])
end
