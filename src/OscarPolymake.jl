module OscarPolymake
using Nemo
using Polymake
#using AbstractAlgebra  #doing Nemo here produces an error
#using Nemo
using LinearAlgebra

export polyhedron, dual_polyhedron

mutable struct HomogeneousPolyhedron #an affine or metric poly object (for me)
  P::Polymake.pm_perl_ObjectAllocated
  tmp::Type
  isbounded::Int #0 -> not known, 1 -> bounded, 2 -> unbounded
  function HomgeneousPolyhedron(t::Type) 
    r = new()
    r.tmp = t
    return r
  end
end

mutable struct Polyhedron #a real polymake polyhedron
  P::HomogeneousPolyhedron
  function Polyhedron(t::Type) 
    r = new()
    r.tmp = t
    return r
  end
end

function show(io::IO, P::Polyhedron)
  println(io, "a (metric) polytope of type $(r.tmp)")
end

# The following is motivated by standard conventions in linear programming.
# ... add ref ...
# CURRENTLY ONLY INTEGER COEFFICIENTS

# P(A,b) = { x : A x <= b }
function polyhedron(A::T, b::T) where {T <: MatElem{<:RingElem}}
  t = typeof(A[1,1])
  #here BigInt, Integer, (fmpz, fmpq) -> Rational
  #     nf_elem quad real field: -> QuadraticExtension
  #     float -> Float
  #     mpfr, BigFloat -> AccurateFloat
  #    
  bA = Array{BigInt, 2}([b -A])
  p = Polymake.perlobj( "Polytope<Rational>", Dict("INEQUALITIES" => bA))
  H = HomogeneousPolyhedron(t)
  H.p = p
  P = Polyhedron(t)
  P.P = H
  return P
end

# P^*(A,c) = { y : y A = c, y >= 0 }
function dual_polyhedron(A::T, c::T) where {T <: MatElem{<:RingElem}}
  t = typeof(A[1,1])
  #here BigInt, Integer, (fmpz, fmpq) -> Rational
  #     nf_elem quad real field: -> QuadraticExtension
  #     float -> Float
  #     mpfr, BigFloat -> AccurateFloat
  #    
  cA = Array{BigInt, 2}([c -transpose(A)])
  m = rows(A)
  nonnegative = [zeros(BigInt,m,1)  Matrix{BigInt}(I,m,m)]
  p_star = Polymake.perlobj( "Polytope<Rational>", Dict("EQUATIONS" => cA, "INEQUALITIES" => nonnegative))
  H = HomogeneousPolyhedron(t)
  H.P = p_star
  P_star = Polyhedron(t)
  P_star.P = H
  return P_star
end

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

greet() = print("Hello World!")

end # module

#=

struct OscarRing
  plus::Ptr{nothing}
  ...
  function OscarRing(R::Ring)
    r = new()
    r.plus = @cfunction(+, elem_type(R), ...)
    r.zero = @cfunction(() -> R(0), ....)
  end
end

pm.init_generic(OscarRing(K))

generic_elem
  ptr::julia
  ptr::struct


=#
