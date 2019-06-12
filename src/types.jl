@doc Markdown.doc"""
    HomogeneousPolyhedron(A)

The homogeneous polyhedron defined by the inequalities $ A x ≥ 0$.
"""
struct HomogeneousPolyhedron #
    P::Polymake.pm_perl_ObjectAllocated
    boundedness::Symbol # Values: :unknown, :bounded, :unbounded
end
function HomogeneousPolyhedron(P::Polymake.pm_perl_ObjectAllocated)
    HomogeneousPolyhedron(P, :unknown)
end
function HomogeneousPolyhedron(bA)
  p = Polymake.perlobj("polytope::Polytope<Rational>", INEQUALITIES=matrix_for_polymake(bA))
  return HomogeneousPolyhedron(p)
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

struct LinearProgram
   polymake_lp::Polymake.pm_perl_ObjectAllocated
   function LinearProgram(P::Polyhedron, objective::Vector{Rational})
      ambDim = Polymake.polytope.AMBIENT_DIM(P.homogeneous_polyhedron.P);
      if size(objective)[1] != ambDim
         error("objective has wrong dimension.")
      end
      lp = Polymake.@pm Polytope.LinearProgram(:LINEAR_OBJECTIVE=>homogenize(o, 0))
      P.LP = lp
      new(lp)
   end
end

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
