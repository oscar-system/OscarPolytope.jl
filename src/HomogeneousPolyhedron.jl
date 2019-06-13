@doc Markdown.doc"""
    HomogeneousPolyhedron(A)

The homogeneous polyhedron defined by the inequalities $ A x ≥ 0$.
"""
struct HomogeneousPolyhedron #
    polymakePolytope::Polymake.pm_perl_ObjectAllocated
    boundedness::Symbol # Values: :unknown, :bounded, :unbounded
end
function HomogeneousPolyhedron(polymakePolytope::Polymake.pm_perl_ObjectAllocated)
    HomogeneousPolyhedron(polymakePolytope, :unknown)
end
function HomogeneousPolyhedron(bA)
  p = Polymake.perlobj("polytope::Polytope<Rational>", INEQUALITIES=matrix_for_polymake(bA))
  return HomogeneousPolyhedron(p)
end

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, H::HomogeneousPolyhedron)
   if(property_is_computed(H, :INEQUALITIES))
      print(io, "Homogeneous polyhedron given by { x | A x ≥ 0 } where \n")
      print(io, "\nA = \n")
      Base.print_array(io, H.polymakePolytope.INEQUALITIES)
   end
end

function Base.isequal(H0::HomogeneousPolyhedron, H1::HomogeneousPolyhedron)
   if(! Polymake.Polytope.included_polyhedra(H0.polymakePolytope, H1.polymakePolytope))
      return false
   end
   if(! Polymake.Polytope.included_polyhedra(H1.polymakePolytope, H0.polymakePolytope))
      return false
   end
   return true
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
function dim(H::HomogeneousPolyhedron)
   return H.polymakePolytope.DIM
end

"""
   vertices(H::HomogeneousPolyhedron)

Returns the vertices of a polyhedron.
"""
function vertices(H::HomogeneousPolyhedron)
   return transpose(H.polymakePolytope.VERTICES)
end

"""
   lineality_space(H::HomogeneousPolyhedron)

Returns a basis of the lineality space of a polyhedron.
"""
function lineality_space(H::HomogeneousPolyhedron)
   return transpose(H.polymakePolytope.LINEALITY_SPACE)
end

"""
   facets(H::HomogeneousPolyhedron)

Returns the facets of a polyhedron.
"""
function facets(H::HomogeneousPolyhedron)
   return H.polymakePolytope.FACETS
end

###############################################################################
###############################################################################
### Standard constructions
###############################################################################
###############################################################################
@doc Markdown.doc"""
   homogeneous_cube(d [, u, l])

Construct the $[-1,1]$-cube in dimension $d$. If $u$ and $l$ are given, the $[l,u]$-cube in dimension $d$ is returned.
"""
function homogeneous_cube(d)
   C = Polymake.Polytope.cube(d)
   return HomogeneousPolyhedron(C)
end
function homogeneous_cube(d, u, l)
   C = Polymake.Polytope.cube(d, u, l)
   return HomogeneousPolyhedron(C)
end
