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
