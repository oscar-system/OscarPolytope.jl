function ==(P0::Polyhedron, P1::Polyhedron)
   return P0.homogeneous_polyhedron == P1.homogeneous_polyhedron
end

@doc Markdown.doc"""
    convex_hull(V [, R [, L]])

The polytope given as the convex hull of the columns of V. Optionally, rays (R)
and generators of the lineality space (L) can be given as well.

see Def. 2.11 and Def. 3.1.
""" function convex_hull(V::AbstractVecOrMat)
   p = polytope.Polytope{Rational}(POINTS = homogenize(transpose(V), 1))
   return Polyhedron(HomogeneousPolyhedron(p))
end
function convex_hull(V::AbstractVecOrMat, R::AbstractVecOrMat)
   p = polytope.Polytope{Rational}(POINTS = vcat(
      homogenize(transpose(V), 1),
      homogenize(transpose(R), 0),
   ))
   return Polyhedron(HomogeneousPolyhedron(p))
end
function convex_hull(
   V::AbstractVecOrMat,
   R::AbstractVecOrMat,
   L::AbstractVecOrMat,
)
   p = polytope.Polytope{Rational}(
      POINTS = vcat(homogenize(transpose(V), 1), homogenize(transpose(R), 0)),
      INPUT_LINEALITY = homogenize(transpose(L), 0),
   )
   return Polyhedron(HomogeneousPolyhedron(p))
end
###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, P::Polyhedron)
   if iscomputed(P, :VERTICES)
      println(
         io,
         "Polyhedron given as the convex hull of the columns of V, where\nV = ",
      )
      Base.print_array(io, vertices(P))
      R = rays(P)
      if size(R, 2) > 0
         println(io, "\n\nwith rays given as the columns of R, where\nR =")
         Base.print_array(io, R)
      end
      L = lineality_space(P)
      if size(L, 2) > 0
         println(
            io,
            "\n\nwith lineality space minimally generated by the columns of L, where\nL =",
         )
         Base.print_array(io, L)
      end
      return
   elseif iscomputed(P, :INEQUALITIES)
      ineq = P.homogeneous_polyhedron.polymakePolytope.INEQUALITIES
      println(io, "Polyhedron given by { x | A x ≤ b } where ")
      println(io, "\nA = ")
      Base.print_array(io, -dehomogenize(ineq))
      println(io, "\n\nb = ")
      Base.print_array(io, ineq[:, 1])
   else
      println(
         io,
         "A Polyhedron with neither vertex nor face representation computed.",
      )
   end
end

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################
"""
   dim(P::Polyhedron)

Returns the dimension of a polyhedron.
"""
polytope.dim(P::Polyhedron) = dim(P.homogeneous_polyhedron)

"""
   ambient_dim(P::Polyhedron)

Returns the ambient dimension of a polyhedron.
"""
polytope.ambient_dim(P::Polyhedron) = ambient_dim(P.homogeneous_polyhedron)

"""
   vertices(P::Polyhedron)

Returns the vertices of a polyhedron in column-major format.
"""
function vertices(P::Polyhedron)
   verts_in_rows = transpose(vertices(P.homogeneous_polyhedron))
   return transpose(dehomogenize(verts_in_rows, vertex_indices(verts_in_rows)))
end

"""
   rays(P::Polyhedron)

Returns minimal set of generators of the cone of unbounded directions of a polyhedron in column-major format.
"""
function rays(P::Polyhedron)
   verts_in_rows = transpose(vertices(P.homogeneous_polyhedron))
   return transpose(dehomogenize(verts_in_rows, ray_indices(verts_in_rows)))
end

"""
   lineality_space(P::Polyhedron)

Returns a basis of the lineality space of a polyhedron in column-major format.
"""
function lineality_space(P::Polyhedron)
   lineality_in_rows = transpose(lineality_space(P.homogeneous_polyhedron))
   return transpose(dehomogenize(
      lineality_in_rows,
      ray_indices(lineality_in_rows),
   ))
end

"""
   facets(P::Polyhedron)

Returns the facets of a polyhedron in column-major format.
"""
facets(P::Polyhedron) =
   decompose_hdata(transpose(facets(P.homogeneous_polyhedron)))

###############################################################################
###############################################################################
### Standard constructions
###############################################################################
###############################################################################
@doc Markdown.doc"""
    dual_polyhedron(A, c)

The (metric) polyhedron defined by

$$P^*(A,c) = \{ y | yA = c, y ≥ 0 \}.$$

see Theorem 4.11.
""" function dual_polyhedron(A, c)
   m, n = size(A)
   cA = transpose([c'; -LinearAlgebra.transpose(A)])
   nonnegative = [zeros(eltype(A), n, 1) LinearAlgebra.I]
   P_star = polytope.Polytope{Rational}(
      EQUATIONS = cA,
      INEQUALITIES = nonnegative,
   )
   H = HomogeneousPolyhedron(P_star)
   return Polyhedron(H)
end

@doc Markdown.doc"""
   cube(d [, u, l])

Construct the $[-1,1]$-cube in dimension $d$. If $u$ and $l$ are given, the $[l,u]$-cube in dimension $d$ is returned.
""" cube(d) = Polyhedron(homogeneous_cube(d))
cube(d, u, l) = Polyhedron(homogeneous_cube(d, u, l))
