@doc Markdown.doc"""
    Polyhedron(A, b)

The (metric) polyhedron defined by

$$P(A,b) = \{ x |  Ax ≤ b \}.$$

see Def. 3.35 and Section 4.1.
"""

struct Polyhedron #a real polymake polyhedron
    homogeneous_polyhedron::HomogeneousPolyhedron
end
Polyhedron(A, b) = Polyhedron(HomogeneousPolyhedron([b -A]))
Polyhedron(pmp::Polymake.pm_perl_ObjectAllocated) = Polyhedron(HomogeneousPolyhedron(pmp))

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

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################
"""
   vertices(P::Polyhedron)

Returns the vertices of a polyhedron.
"""
function vertices(P::Polyhedron)
   result = P.homogeneous_polyhedron.polymakePolytope.VERTICES
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
   result = P.homogeneous_polyhedron.polymakePolytope.VERTICES
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
   result = P.homogeneous_polyhedron.polymakePolytope.LINEALITY_SPACE
   selectedRows = Int[]
   nrows = Polymake.rows(result)
   for i=1:nrows
      if result[i,1] == 0
         push!(selectedRows, i)
      end
   end
   transpose(result[selectedRows, 2:end])
end

###############################################################################
###############################################################################
### Standard constructions
###############################################################################
###############################################################################
@doc Markdown.doc"""
    DualPolyhedron(A, c)

The (metric) polyhedron defined by

$$P^*(A,c) = \{ y | yA = c, y ≥ 0 \}.$$

see Theorem 4.11.
"""
function DualPolyhedron(A, c)
    #here BigInt, Integer, (fmpz, fmpq) -> Rational
    #     nf_elem quad real field: -> QuadraticExtension
    #     float -> Float
    #     mpfr, BigFloat -> AccurateFloat
    #
    m, n = size(A)
    cA = matrix_for_polymake([c -LinearAlgebra.transpose(A)])
    nonnegative = [zeros(BigInt, m, 1)  LinearAlgebra.I]
    P_star = Polymake.perlobj("polytope::Polytope<Rational>", EQUATIONS=cA, INEQUALITIES=nonnegative)
    H = HomogeneousPolyhedron(P_star)
    return Polyhedron(H)
end

@doc Markdown.doc"""
   cube(d [, u, l])

Construct the $[-1,1]$-cube in dimension $d$. If $u$ and $l$ are given, the $[l,u]$-cube in dimension $d$ is returned.
"""
function cube(d)
   C = Polymake.Polytope.cube(d)
   return Polyhedron(C)
end
function cube(d, u, l)
   C = Polymake.Polytope.cube(d, u, l)
   return Polyhedron(C)
end
