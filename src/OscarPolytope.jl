module OscarPolytope

import LinearAlgebra, Markdown, Nemo, Polymake

export Polyhedron, HomogenousPolyhedron

# export polyhedron, dual_polyhedron, homogenous_polyhedron,
       # boundary_points_matrix, inner_points_matrix

#to add the (common) interface... still should be renamed to nrows
# Nemo.rows(A::Polymake.pm_MatrixAllocated) = Int(size(A)[1])
# Nemo.cols(A::Polymake.pm_MatrixAllocated) = Int(size(A)[2])


struct HomogenousPolyhedron #an affine or metric poly object (for me)
    P::Polymake.pm_perl_ObjectAllocated
    boundedness::Symbol # Values: :unknown, :bounded, :unbounded
end
function HomogenousPolyhedron(P::Polymake.pm_perl_ObjectAllocated)
    HomogenousPolyhedron(P, :unknown)
end
@doc Markdown.doc"""
    HomogenousPolyhedron(A)

# The homogenous polyhedron defined by the inequalities $ A x ≥ 0$.
"""
function HomogenousPolyhedron(bA::AbstractMatrix)
  #here BigInt, Integer, (fmpz, fmpq) -> Rational
  #     nf_elem quad real field: -> QuadraticExtension
  #     float -> Float
  #     mpfr, BigFloat -> AccurateFloat
  #
  p = Polymake.perlobj("Polytope<Rational>", Dict("INEQUALITIES" => Rational{BigInt}.(bA)))
  return HomogenousPolyhedron(p)
end

function Base.show(io::IO, H::HomogenousPolyhedron)
    print(io, "Homogeneous polyhedron given by { x | A x ≤ 0 } where \n")
    print(io, "\nA = \n")
    Base.print_array(io, H.P.INEQUALITIES)
end

@doc Markdown.doc"""
    Polyhedron(A, b)

The (metric) polyhedron defined by

$$P(A,b) = \{ x \, | \, Ax \le b \}.$$
"""
struct Polyhedron #a real polymake polyhedron
    homogenous_polyhedron::HomogenousPolyhedron
end
Polyhedron(A, b) = Polyhedron(HomogenousPolyhedron([b -A]))

function Base.show(io::IO, P::Polyhedron)
    ineq = P.homogenous_polyhedron.P.INEQUALITIES
    print(io, "Polyhedron given by { x | A x ≤ b } where \n")
    print(io, "\nA = \n")
    Base.print_array(io, ineq[:,2:end])
    print(io, "\n\nb = \n")
    Base.print_array(io, ineq[:,1])
end


# function homogenous_polyhedron(bA::T) where {T <: MatElem{<:RingElem}}
#     return homogenous_polyhedron(Matrix{BigInt}(bA))
# end
#
#
# # P^*(A,c) = { y : y A = c, y >= 0 }
# @doc Markdown.doc"""
#     dual_polyhedron(A, c) -> Polyhedron
#
# > The dual (metric) polyhedron defined by $yA = c$ and $ y \ge 0$.
# """
# function dual_polyhedron(A::T, c::T) where {T <: MatElem{<:RingElem}}
#   t = typeof(A[1,1])
#   #here BigInt, Integer, (fmpz, fmpq) -> Rational
#   #     nf_elem quad real field: -> QuadraticExtension
#   #     float -> Float
#   #     mpfr, BigFloat -> AccurateFloat
#   #
#   cA = Array{BigInt, 2}([c -transpose(A)])
#   m = rows(A)
#   nonnegative = [zeros(BigInt,m,1)  Matrix{BigInt}(I,m,m)]
#   p_star = Polymake.perlobj( "Polytope<Rational>", Dict("EQUATIONS" => cA, "INEQUALITIES" => nonnegative))
#   H = HomogenousPolyhedron(t)
#   H.P = p_star
#   P_star = Polyhedron(H)
#   return P_star
# end
#
# _polytope(P::Polyhedron) = P.P.P
# _polytope(H::HomogenousPolyhedron) = H.P
#
# #we don't have points yet, so I can only return the matrix.
# # the polyhedron is not homogenous, so I strip the 1st entry
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
