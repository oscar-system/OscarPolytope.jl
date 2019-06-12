using Test
using OscarPolytope

@testset "convex_hull" begin
pts = [1 0 0; 0 0 1]
Q = convex_hull(pts)
V = vertices(Q)
@test V == pts 
end
