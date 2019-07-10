function HomogeneousPolyhedron(polymakePolytope::Polymake.pm_perl_ObjectAllocated)
    HomogeneousPolyhedron(polymakePolytope, :unknown)
end
function HomogeneousPolyhedron(bA)
  p = Polymake.@pm Polytope.Polytope{Rational}(:INEQUALITIES=>matrix_for_polymake(bA))
  return HomogeneousPolyhedron(p)
end

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, H::HomogeneousPolyhedron)
   if property_is_computed(H, :INEQUALITIES)
      print(io, "Homogeneous polyhedron given by { x | A x ≥ 0 } where \n")
      print(io, "\nA = \n")
      Base.print_array(io, H.polymakePolytope.INEQUALITIES)
   else
      println(io, "Homogeneous polyhedron defined as convex hull of vertices and rays.")
      println(io, "The hyperplane description has not been computed yet.")
   end
end

function ==(H0::HomogeneousPolyhedron, H1::HomogeneousPolyhedron)
   return Polytope.included_polyhedra(H0.polymakePolytope, H1.polymakePolytope) &&
      Polytope.included_polyhedra(H1.polymakePolytope, H0.polymakePolytope)
end

###############################################################################
###############################################################################
### Access properties
###############################################################################
###############################################################################
"""
   dim(H::HomogeneousPolyhedron)

Returns the dimension of a polyhedron.
"""
Polytope.dim(H::HomogeneousPolyhedron) = Polytope.dim(H.polymakePolytope)

"""
   ambient_dim(H::HomogeneousPolyhedron)

Returns the ambient dimension of a polyhedron.
"""
Polytope.ambient_dim(H::HomogeneousPolyhedron) = Polytope.ambient_dim(H.polymakePolytope)

"""
   vertices(H::HomogeneousPolyhedron)

Returns the vertices of a polyhedron.
"""
vertices(H::HomogeneousPolyhedron) = transpose(H.polymakePolytope.VERTICES)

"""
   lineality_space(H::HomogeneousPolyhedron)

Returns a basis of the lineality space of a polyhedron.
"""
lineality_space(H::HomogeneousPolyhedron) = transpose(H.polymakePolytope.LINEALITY_SPACE)

"""
   facets(H::HomogeneousPolyhedron)

Returns the facets of a polyhedron.
"""
facets(H::HomogeneousPolyhedron) = H.polymakePolytope.FACETS

###############################################################################
###############################################################################
### Standard constructions
###############################################################################
###############################################################################
@doc Markdown.doc"""
   homogeneous_cube(d [, u, l])

Construct the $[-1,1]$-cube in dimension $d$. If $u$ and $l$ are given, the $[l,u]$-cube in dimension $d$ is returned.
"""
homogeneous_cube(d) = HomogeneousPolyhedron(Polytope.cube(d))
homogeneous_cube(d, u, l) = HomogeneousPolyhedron(Polytope.cube(d, u, l))
