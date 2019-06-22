using Test
using OscarPolytope

pts = [1 0 0; 0 0 1]
Q0 = convex_hull(pts)
Q1 = convex_hull(pts, [1;1]);
Q2 = convex_hull(pts, [1;1], [1;1]);
C0 = cube(2)
C1 = cube(2,1,0)

@testset "OscarPolytope" begin

@testset "(de)homogenize/augment" begin
   dehomogenize, homogenize = OscarPolytope.dehomogenize, OscarPolytope.homogenize
   pm = OscarPolytope.Polymake
   m = [1 2 3; 4 5 6]
   @test dehomogenize(homogenize(m, 0//1)) == m
   @test dehomogenize(homogenize(m)) == m
   @test dehomogenize(homogenize(pm.pm_Matrix(m))) == m
   @test dehomogenize(homogenize(pm.pm_Matrix{pm.pm_Integer}(m))) isa pm.pm_Matrix{pm.pm_Integer}
   @test dehomogenize(homogenize(pm.pm_Matrix{pm.pm_Rational}(m))) isa pm.pm_Matrix{pm.pm_Rational}

   v = [1,2,3]
   @test dehomogenize(homogenize(v, 1//1)) == v
   @test dehomogenize(homogenize(v)) == v
   @test dehomogenize(homogenize(pm.pm_Vector{pm.pm_Integer}(v))) == v
   @test dehomogenize(homogenize(pm.pm_Vector{pm.pm_Integer}(v))) isa pm.pm_Vector{pm.pm_Integer}
   @test dehomogenize(homogenize(pm.pm_Vector{pm.pm_Rational}(v))) isa pm.pm_Vector{pm.pm_Rational}

   augment = OscarPolytope.augment

   @test augment(m, [9,10]) == [9 1 2 3; 10 4 5 6]
   @test augment(pm.pm_Matrix(m), [9,10]) == [9 1 2 3; 10 4 5 6]
end

@testset "conformance tests" begin
   @test typeof(Q0) == Polyhedron
   @test typeof(Q1) == Polyhedron
   @test typeof(Q2) == Polyhedron
   @test typeof(C0) == Polyhedron
   @test typeof(C1) == Polyhedron
   @test typeof(Q0 == Q0) == Bool
   @test typeof(Q0 == Q1) == Bool
   @test typeof(Q0 != Q0) == Bool
   @test typeof(Q0 != Q1) == Bool
   @test Q0 != Q1
   @test C0 != C1
   @test C0 == C0
   @test typeof(dim(Q0)) == Int
   @test typeof(ambient_dim(Q0)) == Int
   @test Q2 == convex_hull(vertices(Q2), rays(Q2), lineality_space(Q2))
end

@testset "convex_hull" begin
   @test size(vertices(Q0)) == (2,3)
   @test size(vertices(Q1)) == (2,3)
   @test size(rays(Q1)) == (2,1)
   @test size(lineality_space(Q1)) == (2,0)
   @test size(vertices(Q2)) == (2,2)
   @test size(rays(Q2)) == (2,0)
   @test size(lineality_space(Q2)) == (2,1)
   @test dim(Q0) == 2
   @test dim(Q1) == 2
   @test dim(Q2) == 2
   @test ambient_dim(Q0) == 2
   @test ambient_dim(Q1) == 2
   @test ambient_dim(Q2) == 2
end

@testset "standard constructions" begin
   @test size(vertices(C0)) == (2,4)
   @test C0 == convex_hull(vertices(C0))
end

@testset "LinearProgram" begin
   A = [-1 0; 0 -1; 1 0; 0 1]
   b = [0; 0; 1; 1]
   objective = [1;1]
   primal = PrimalProgram(objective,A,b)
   dual = DualProgram(objective,A,b)
   @test typeof(primal) == LinearProgram
   @test minimal_value(primal) == 0
   @test maximal_value(primal) == 2
   @test minimal_vertex(primal) == [0;0]
   @test maximal_vertex(primal) == [1;1]
   @test maximal_value(primal) == minimal_value(dual)
   @test minimal_vertex(dual) == [0,0,1,1]
end

end # of @testset "OscarPolytope"
