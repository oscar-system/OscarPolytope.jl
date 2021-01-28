struct LinearProgram
   feasible_region::Polyhedron
   polymake_lp::Polymake.BigObjectAllocated
   function LinearProgram(P::Polyhedron, objective::AbstractVector)
      ambDim = ambient_dim(P)
      size(objective, 1) == ambDim || error("objective has wrong dimension.")
      lp = Polymake.polytope.LinearProgram(LINEAR_OBJECTIVE=>homogenize(objective, 0))
      pm_polytope(P).LP = lp
      new(P, lp)
   end
end



"""
   PrimalProgram(c, A, b)

Constructs the primal linear program max{cx | Ax<= b}.

see Def. 4.10
"""
PrimalProgram(c, A, b) = LinearProgram(Polyhedron(A,b), c)

# """
#    DualProgram(c, A, b)
#
# Constructs the dual linear program min{yb | yA=c, y>=0}.
#
# see Def. 4.10
# """
# DualProgram(c, A, b) = LinearProgram(DualPolyhedron(A,c), b)

minimal_vertex(lp::LinearProgram) = dehomogenize(lp.polymake_lp.MINIMAL_VERTEX)
maximal_vertex(lp::LinearProgram) = dehomogenize(lp.polymake_lp.MAXIMAL_VERTEX)

minimal_value(lp::LinearProgram) = lp.polymake_lp.MINIMAL_VALUE
maximal_value(lp::LinearProgram) = lp.polymake_lp.MAXIMAL_VALUE
