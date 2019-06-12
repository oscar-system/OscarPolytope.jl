module OscarPolytope

import LinearAlgebra, Markdown, Nemo, Polymake

export Polyhedron, DualPolyhedron, HomogeneousPolyhedron, vertices, rays, LinearProgram, minimal_vertex, minimal_value, maximal_vertex, maximal_value, convex_hull, property_is_computed, lineality_space

include("types.jl")

# export polyhedron, dual_polyhedron, homogeneous_polyhedron,
       # boundary_points_matrix, inner_points_matrix
matrix_for_polymake(x::Nemo.fmpz_mat) = Matrix{BigInt}(x)
matrix_for_polymake(x::Nemo.fmpq_mat) = Matrix{Rational{BigInt}}(x)
matrix_for_polymake(x::AbstractMatrix{<:Integer}) = x
matrix_for_polymake(x::AbstractMatrix{<:Rational{<:Integer}}) = x

#here BigInt, Integer, (fmpz, fmpq) -> Rational
#     nf_elem quad real field: -> QuadraticExtension
#     float -> Float
#     mpfr, BigFloat -> AccurateFloat

function Base.show(io::IO, H::HomogeneousPolyhedron)
    print(io, "Homogeneous polyhedron given by { x | A x ≥ 0 } where \n")
    print(io, "\nA = \n")
    Base.print_array(io, H.P.INEQUALITIES)
end

@doc Markdown.doc"""
    convex_hull(V [, R [, L]])

The polytope given as the convex hull of the columns of V. Optionally, rays (R)
and generators of the lineality space (L) can be given as well.

see Def. 2.11 and Def. 3.1.
"""
function convex_hull(V)
   p = Polymake.perlobj("polytope::Polytope<Rational>", POINTS=homogenize(transpose(V), 1))
   return Polyhedron(HomogeneousPolyhedron(p))
end
function convex_hull(V, R)
   p = Polymake.perlobj("polytope::Polytope<Rational>", POINTS=vcat(homogenize(transpose(V), 1), homogenize(transpose(R), 0)))
   return Polyhedron(HomogeneousPolyhedron(p))
end
function convex_hull(V, R, L)
   p = Polymake.perlobj("polytope::Polytope<Rational>", POINTS=vcat(homogenize(transpose(V), 1), homogenize(transpose(R), 0)), INPUT_LINEALITY=homogenize(transpose(L),0))
   return Polyhedron(HomogeneousPolyhedron(p))
end

augment(vec::AbstractVector{T}, val) where T = vcat(T(val), vec)
augment(mat::AbstractMatrix, vec::AbstractVector) = hcat(vec, mat)

homogenize(::Type{T}, vec::AbstractVector, val::T=T(0)) where T = augment(T.(vec), val)
homogenize(::Type{T}, mat::AbstractMatrix, val::T=T(1)) where T = augment(T.(mat), repeat([val], size(mat,1)))
homogenize(mat::AbstractVecOrMat{T}, val=T(1)) where T = homogenize(T, mat,
val)

dehomogenize(vec::AbstractVector) = vec[2:end]
dehomogenize(mat::AbstractMatrix) = mat[:, 2:end]

function minimal_vertex(lp::LinearProgram)
   result = lp.polymake_lp.MINIMAL_VERTEX
   result[2:end]
end
function minimal_value(lp::LinearProgram)
   lp.polymake_lp.MINIMAL_VALUE
end

function maximal_vertex(lp::LinearProgram)
   result = lp.polymake_lp.MAXIMAL_VERTEX
   result[2:end]
end
function maximal_value(lp::LinearProgram)
   lp.polymake_lp.MAXIMAL_VALUE
end

function property_is_computed(P::Polymake.pm_perl_ObjectAllocated, S::Symbol)
   pv = Polymake.internal_call_method("lookup", P, Any[string(S)])
   return nothing != Polymake.convert_from_property_value(pv)
end
function property_is_computed(P::Polyhedron, S::Symbol)
   return property_is_computed(P.homogeneous_polyhedron.P, S)
end

function Base.show(io::IO, P::Polyhedron)
   if(property_is_computed(P, :INEQUALITIES))
      ineq = P.homogeneous_polyhedron.P.INEQUALITIES
      print(io, "Polyhedron given by { x | A x ≤ b } where \n")
      print(io, "\nA = \n")
      Base.print_array(io, -ineq[:,2:end])
      print(io, "\n\nb = \n")
      Base.print_array(io, ineq[:,1])
   end
end


"""
   vertices(P::Polyhedron)

Returns the vertices of a polyhedron.
"""
function vertices(P::Polyhedron)
   result = P.homogeneous_polyhedron.P.VERTICES
   selectedRows = Int[]
   nrows = Polymake.rows(result)
   for i=1:nrows
      if result[i,1] == 1
         push!(selectedRows, i)
      end
   end
   transpose(result[selectedRows, 2:end])
end

"""
   rays(P::Polyhedron)

Returns the generators of the cone of unbounded directions of a polyhedron.
"""
function rays(P::Polyhedron)
   result = P.homogeneous_polyhedron.P.VERTICES
   selectedRows = Int[]
   nrows = Polymake.rows(result)
   for i=1:nrows
      if result[i,1] == 0
         push!(selectedRows, i)
      end
   end
   transpose(result[selectedRows, 2:end])
end

"""
   lineality_space(P::Polyhedron)

Returns the generators of the lineality space of a polyhedron.
"""
function lineality_space(P::Polyhedron)
   result = P.homogeneous_polyhedron.P.LINEALITY_SPACE
   selectedRows = Int[]
   nrows = Polymake.rows(result)
   for i=1:nrows
      if result[i,1] == 0
         push!(selectedRows, i)
      end
   end
   transpose(result[selectedRows, 2:end])
end

#
# #we don't have points yet, so I can only return the matrix.
# # the polyhedron is not homogeneous, so I strip the 1st entry
# #TODO: since P is given by Ax <= b, should the "points" be rows or columns?
#
# @doc Markdown.doc"""
#     boundary_points_matrix(P::Polyhedron) -> fmpz_mat
# > Given a bounded polyhedron, return the coordinates of all integer points
# > its boundary as rows in a matrix.
# """
# function boundary_points_matrix(P::Polyhedron)
#   p = _polytope(P)
#   out = Polymake.give(p, "BOUNDARY_LATTICE_POINTS")
#   res = zero_matrix(FlintZZ, rows(out), cols(out)-1)
#   for i=1:rows(out)
#     @assert out[i,1] == 1
#     for j=1:cols(out)-1
#       res[i,j] = out[i, j+1]
#     end
#   end
#   return res
# end
#
# @doc Markdown.doc"""
#     inner_points_matrix(P::Polyhedron) -> fmpz_mat
# > Given a bounded polyhedron, return the coordinates of all integer points
# > in its interior as rows in a matrix.
# """
# function inner_points_matrix(P::Polyhedron)
#   p = _polytope(P)
#   out = Polymake.give(p, "INTERIOR_LATTICE_POINTS")
#   res = zero_matrix(FlintZZ, rows(out), cols(out)-1)
#   for i=1:rows(out)
#     @assert out[i,1] == 1
#     for j=1:cols(out)-1
#       res[i,j] = out[i, j+1]
#     end
#   end
#   return res
# end

#=
function polyhedron(V::Array{Point, 1}, C::Array{Point, 1}=Array{Point, 1}(UndefInitializer(), 0), L::Array{Point, 1} = Array{Point, 1}(UndefInitializer(), 0))# only for Homogeneous version; Other::Dict{String, Any} = Dict{String, Any}())
  d = Dict{String, Any}("INEQUALITIES" => bA)
#  for (k,v) = Other
#    d[k] = v
#  end

#  POINTS = [ 1 V; 0 C], INPUT_LINEALITIES = L
  p = Polymake.perlobj("Polytope<Rational>", d)
  P = Polyhedron(t)
  P.P = p
  return P
end

function polyhedron(V::MatElem, C::MatElem, L::MatElem)
end

function convex_hull(P::polyhedron, Q::Polyhedron)
end
function convex_hull(P::polyhedron, Q::Point)
end
+(P, Q)
intersect(P, Q)
product
join
free_sum

lattice_points
isbounded
isunbounded
lattice_points_generators
volume
isempty


=#
end # module
