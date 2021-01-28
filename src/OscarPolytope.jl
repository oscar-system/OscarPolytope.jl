module OscarPolytope

import LinearAlgebra, Markdown, Polymake
import Base: ==
# import Polymake.polytope: dim, ambient_dim
const Polytope = Polymake.polytope

export Polyhedron
export Polyhedron,
    dim,
    ambient_dim,
    vertices,
    facets,
    rays,
    lineality_space,
    cube,
    convex_hull,
    pm_polytope

include("helpers.jl")
include("Polyhedron.jl")

end # module

###############
## OLD STUFF ##
###############

# include("LinearProgram.jl")

# export polyhedron, dual_polyhedron, homogeneous_polyhedron,
# boundary_points_matrix, inner_points_matrix

#here BigInt, Integer, (fmpz, fmpq) -> Rational
#     nf_elem quad real field: -> QuadraticExtension
#     float -> Float
#     mpfr, BigFloat -> AccurateFloat




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
