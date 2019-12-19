"""
   PrimalProgram(c, A, b)

Constructs the primal linear program max{cx | Ax<= b}.

see Def. 4.10
"""
PrimalProgram(c, A, b) = polytope.LinearProgram(Polyhedron(A,b), c)

"""
   DualProgram(c, A, b)

Constructs the dual linear program min{yb | yA=c, y>=0}.

see Def. 4.10
"""
DualProgram(c, A, b) = polytope.LinearProgram(DualPolyhedron(A,c), b)

minimal_vertex(lp::LinearProgram) = dehomogenize(lp.polymake_lp.MINIMAL_VERTEX)
maximal_vertex(lp::LinearProgram) = dehomogenize(lp.polymake_lp.MAXIMAL_VERTEX)

minimal_value(lp::LinearProgram) = lp.polymake_lp.MINIMAL_VALUE
maximal_value(lp::LinearProgram) = lp.polymake_lp.MAXIMAL_VALUE
