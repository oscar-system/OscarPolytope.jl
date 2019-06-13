struct LinearProgram
   feasible_region::Polyhedron
   polymake_lp::Polymake.pm_perl_ObjectAllocated
   function LinearProgram(P::Polyhedron, objective::AbstractVector)
      ambDim = ambient_dim(P)
      if size(objective)[1] != ambDim
         error("objective has wrong dimension.")
      end
      lp = Polymake.@pm Polytope.LinearProgram(:LINEAR_OBJECTIVE=>homogenize(objective, 0))
      P.homogeneous_polyhedron.polymakePolytope.LP = lp
      new(P, lp)
   end
end

"""
   PrimalProgram(c, A, b)

Constructs the primal linear program max{cx | Ax<= b}.

see Def. 4.10
"""
function PrimalProgram(c, A, b)
   fr = Polyhedron(A,b)
   return LinearProgram(fr, c)
end

"""
   DualProgram(c, A, b)

Constructs the dual linear program min{yb | yA=c, y>=0}.

see Def. 4.10
"""
function DualProgram(c, A, b)
   fr = DualPolyhedron(A,c)
   return LinearProgram(fr, b)
end


function minimal_vertex(lp::LinearProgram)
   result = lp.polymake_lp.MINIMAL_VERTEX
   result[2:end]
end
function minimal_value(lp::LinearProgram)
   lp.polymake_lp.MINIMAL_VALUE
end

function maximal_vertex(lp::LinearProgram)
   result = lp.polymake_lp.MAXIMAL_VERTEX
   result[2:end]
end
function maximal_value(lp::LinearProgram)
   lp.polymake_lp.MAXIMAL_VALUE
end

