@doc Markdown.doc"""
    HomogeneousPolyhedron(A)

The homogeneous polyhedron defined by the inequalities $ A x ≥ 0$.
"""
struct HomogeneousPolyhedron #
    polymakePolytope::Polymake.pm_perl_ObjectAllocated
    boundedness::Symbol # Values: :unknown, :bounded, :unbounded
end


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

struct LinearProgram
   feasible_region::Polyhedron
   polymake_lp::Polymake.pm_perl_ObjectAllocated
   function LinearProgram(P::Polyhedron, objective::AbstractVector)
      ambDim = ambient_dim(P)
      size(objective, 1) == ambDim || error("objective has wrong dimension.")
      lp = Polymake.@pm Polytope.LinearProgram(:LINEAR_OBJECTIVE=>homogenize(objective, 0))
      P.homogeneous_polyhedron.polymakePolytope.LP = lp
      new(P, lp)
   end
end
