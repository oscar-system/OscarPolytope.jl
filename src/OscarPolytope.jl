module OscarPolytope

import LinearAlgebra, Markdown, Nemo, Polymake
import Base: ==

const polytope = Polymake.polytope

import Polymake.polytope: ambient_dim, dim

export Polyhedron,
       HomogeneousPolyhedron,
       LinearProgram,
       ambient_dim,
       augment,
       convex_hull,
       dehomogenize,
       dim,
       dual_program,
       dual_polyhedron,
       facets,
       homogenize,
       lineality_space,
       maximal_vertex,
       maximal_value,
       minimal_vertex,
       minimal_value,
       primal_program,
       rays,
       vertices

include("types.jl")
include("HomogeneousPolyhedron.jl")
include("Polyhedron.jl")
include("helpers.jl")
include("LinearProgram.jl")

end # module
