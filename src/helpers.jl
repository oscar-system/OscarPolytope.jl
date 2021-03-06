matrix_for_polymake(x::Union{Oscar.fmpz_mat,AbstractMatrix{Oscar.fmpz}}) = Matrix{BigInt}(x)
matrix_for_polymake(x::Union{Oscar.fmpq_mat,AbstractMatrix{Oscar.fmpq}}) =
    Matrix{Rational{BigInt}}(x)
matrix_for_polymake(x) = x

function remove_zero_rows(A::Union{Oscar.MatElem,AbstractMatrix})
    A[findall(x->!iszero(x),collect(eachrow(A))),:]
end

function augment(vec::AbstractVector, val)
    s = size(vec)
    res = similar(vec, (s[1] + 1,))
    res[1] = val
    res[2:end] = vec
    return res
end

function augment(mat::AbstractMatrix, vec::AbstractVector)
    s = size(mat)
    res = similar(mat, (s[1], s[2] + 1))
    res[:, 1] = vec
    res[:, 2:end] = mat
    return res
end

homogenize(vec::AbstractVector, val::Number = 0) = augment(vec, val)
homogenize(mat::AbstractMatrix, val::Number = 1) = augment(mat, fill(val, size(mat, 1)))

dehomogenize(vec::AbstractVector) = vec[2:end]
dehomogenize(mat::AbstractMatrix) = mat[:, 2:end]

"""
    stack(A::AbstractVecOrMat, B::AbstractVecOrMat)

Stacks `A` and `B` vertically. The difference to `vcat`is that `AbstractVector`s are always
interpreted as row vectors. Empty vectors are ignored.

## Examples

```
julia> stack([1, 2], [0, 0])
2×2 Array{Int64,2}:
 1  2
 0  0

julia> stack([1 2], [0 0])
2×2 Array{Int64,2}:
 1  2
 0  0

julia> stack([1 2], [0, 0])
2×2 Array{Int64,2}:
 1  2
 0  0

julia> stack([1, 2], [0 0])
2×2 Array{Int64,2}:
 1  2
 0  0

julia> stack([1 2], [])
1×2 Array{Int64,2}:
 1  2
```
"""
stack(A::AbstractMatrix, B::AbstractMatrix) = [A; B]
stack(A::AbstractMatrix, B::AbstractVector) = isempty(B) ? A :  [A; B']
stack(A::AbstractVector, B::AbstractMatrix) = isempty(A) ? B : [A'; B]
stack(A::AbstractVector, B::AbstractVector) = isempty(A) ? B : [A'; B']
#=
function stack(A::Array{Polymake.VectorAllocated{Polymake.Rational},1})
    if length(A)==2
        return(stack(A[1],A[2]))
    end
    M=stack(A[1],A[2])
    for i in 3:length(A)
        M=stack(M,A[i])
    end
    return(M)
end
=#

function property_is_computed(P::Polymake.BigObjectAllocated, S::Symbol)
    pv = Polymake.internal_call_method("lookup", P, Any[string(S)])
    return nothing != Polymake.convert_from_property_value(pv)
end

"""
   decompose_vdata(A::AbstractMatrix)

Given a (homogeneous) polymake matrix split into vertices and rays and dehomogenize.
"""
function decompose_vdata(A::AbstractMatrix)
    vertex_indices = findall(!iszero, view(A, :, 1))
    ray_indices = findall(iszero, view(A, :, 1))
    return (vertices = A[vertex_indices, 2:end], rays = A[ray_indices, 2:end])
end

function decompose_hdata(A)
    (A = -A[:, 2:end], b = A[:, 1])
end
