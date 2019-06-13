struct LinearProgram
   polymake_lp::Polymake.pm_perl_ObjectAllocated
   function LinearProgram(P::Polyhedron, objective::Vector{Rational})
      ambDim = Polymake.polytope.AMBIENT_DIM(P.homogeneous_polyhedron.polymakePolytope);
      if size(objective)[1] != ambDim
         error("objective has wrong dimension.")
      end
      lp = Polymake.@pm Polytope.LinearProgram(:LINEAR_OBJECTIVE=>homogenize(o, 0))
      P.LP = lp
      new(lp)
   end
end

