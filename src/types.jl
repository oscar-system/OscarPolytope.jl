@doc Markdown.doc"""
    HomogeneousPolyhedron

The polyhedron in homogeneous coordinates defined by the inequalities $ A x ≥ b$.
The existing constructors are:
 * `HomogeneousPolyhedron(P::pm_perl_Object)`: constructing a polyhedron from a `Polymake` polytope
 * `HomogeneousPolyhedron(bA::AbstractMatrix)`: constructing a polyhedron from matrix of homogeneous inequalities, i.e. `[b A]`, where matrix `A` describes the inequalities and vector `b` contains their intercepts.
 NOTE: This constructor needs ROW-MAJOR (e.g. polymake) format!
 * `HomogeneousPolyhedron(A::AbstractMatrix, b::AbstractVector)`: constructing a polyhedron from matrix of homogeneous inequalities, i.e. `[b -A]`, where matrix `A` describes the inequalities (column-major) and vector `b` contains their intercepts.
""" struct HomogeneousPolyhedron #
    polymakePolytope::Polymake.pm_perl_ObjectAllocated
    boundedness::Symbol # Values: :unknown, :bounded, :unbounded

    HomogeneousPolyhedron(P::Polymake.pm_perl_Object) = new(P)

    function HomogeneousPolyhedron(
        bA::Union{AbstractMatrix,Nemo.MatrixElem},
        boundedness::Symbol = :unknown,
    )
        bounded_attrs = (:unknown, :bounded, :unbounded)
        @assert boundedness in bounded_attrs "`bounded` attribute may be either of $bounded_attrs"
        return new(polytope.Polytope(:INEQUALITIES => bA), boundedness)
    end

    HomogeneousPolyhedron(A::AbstractMatrix, b::AbstractVector) =
        HomogeneousPolyhedron(transpose([b'; -A]))
end


@doc Markdown.doc"""
    Polyhedron(A, b)

The (metric) polyhedron defined by

$$P(A,b) = \{ x |  Ax ≤ b \}.$$

see Def. 3.35 and Section 4.1. Matrix `A` must be given in column-major format.
""" struct Polyhedron #a real polymake polyhedron
    homogeneous_polyhedron::HomogeneousPolyhedron

    Polyhedron(hp::HomogeneousPolyhedron) = new(hp)

    Polyhedron(pmp::Polymake.pm_perl_Object) = new(HomogeneousPolyhedron(pmp))

    Polyhedron(A, b) = new(HomogeneousPolyhedron(transpose([b'; -A])))

    Polyhedron(bA) = new(HomogeneousPolyhedron(bA))
end

struct LinearProgram
    feasible_region::Polyhedron
    polymake_lp::Polymake.pm_perl_ObjectAllocated
    function LinearProgram(P::Polyhedron, objective::AbstractVector)
        ambDim = ambient_dim(P)
        size(objective, 1) == ambDim || Throw(ArgumentError("objective has wrong dimension."))
        lp = polytope.LinearProgram(:LINEAR_OBJECTIVE => homogenize(
            objective,
            0,
        ))
        P.homogeneous_polyhedron.polymakePolytope.LP = lp
        return new(P, lp)
    end
end
