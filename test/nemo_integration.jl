using OscarPolytope.Polymake
using Nemo

@testset "Nemo Integration" begin
   m = matrix(ZZ, [1 2 3; 4 5 6])
   @test polytope.Polytope(INEQUALITIES=m) isa Polymake.pm_perl_Object

   # @test Polyhedron(m) isa Polyhedron
end
