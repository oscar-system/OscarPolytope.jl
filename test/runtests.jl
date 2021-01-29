using Test
using OscarPolytope
const pm = OscarPolytope.Polymake


@testset "OscarPolytope" begin

    pts = [1 0 0; 0 0 1]'
    Q0 = convex_hull(pts)
    Q1 = convex_hull(pts, [1 1])
    Q2 = convex_hull(pts, [1 1], [1 1])
    C0 = cube(2)
    C1 = cube(2, 1, 0)

    @testset "(de)homogenize" begin
        dehomogenize, homogenize = OscarPolytope.dehomogenize, OscarPolytope.homogenize

        m = [1 2 3; 4 5 6]'
        @test dehomogenize(homogenize(m, 0 // 1)) == m
        @test dehomogenize(homogenize(m)) == m
        @test dehomogenize(homogenize(pm.Matrix(m))) == m
        @test dehomogenize(homogenize(pm.Matrix{pm.Integer}(m))) isa
              pm.Matrix{pm.Integer}
        @test dehomogenize(homogenize(pm.Matrix{pm.Rational}(m))) isa
              pm.Matrix{pm.Rational}
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
        @test vertices(Q0) == vertices(convex_hull(vertices(Q0)))
    end

    @testset "convex_hull" begin
        @test size(vertices(Q0)) == (3, 2)
        @test size(vertices(Q1)) == (3, 2)
        @test size(rays(Q1)) == (1, 2)
        @test size(lineality_space(Q1)) == (0, 2)
        @test size(vertices(Q2)) == (2, 2)
        @test size(rays(Q2)) == (0, 2)
        @test size(lineality_space(Q2)) == (1, 2)
        @test dim(Q0) == 2
        @test dim(Q1) == 2
        @test dim(Q2) == 2
        @test ambient_dim(Q0) == 2
        @test ambient_dim(Q1) == 2
        @test ambient_dim(Q2) == 2
    end

    @testset "standard constructions" begin
        @test size(vertices(C0)) == (4, 2)
        @test_broken C0 == convex_hull(vertices(C0))
    end
end # of @testset "OscarPolytope"
