using Test
using OscarPolytope

@testset "convex_hull" begin
pts = [1 0 0; 0 0 1]
Q0 = convex_hull(pts)
V = vertices(Q0)
@test V == pts 
Q1 = convex_hull(pts, [1;1]);
@test size(vertices(Q1)) == (2,3)
@test size(rays(Q1)) == (2,1)
@test size(lineality_space(Q1)) == (2,0)
Q2 = convex_hull(pts, [1;1], [1;1]);
@test size(vertices(Q2)) == (2,2)
@test size(rays(Q2)) == (2,0)
@test size(lineality_space(Q2)) == (2,1)
end
