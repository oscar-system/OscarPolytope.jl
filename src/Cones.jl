@doc Markdown.doc"""
    Cone(Rays)

     A polyhedral cone, not necessarily pointed, defined by the positive hull
     of the rays, Rays.
""" struct Cone #a real polymake polyhedron
    pm_cone::Polymake.BigObjectAllocated
end
function Cone(Rays::Union{Oscar.MatElem,AbstractMatrix})
    Cone(Polymake.polytope.Cone{Rational}(
        RAYS = matrix_for_polymake(Rays),
    ))
end

"""
    positive_hull(Rays::Union{Oscar.MatElem,AbstractMatrix})

     A polyhedral cone, not necessarily pointed, defined by the positive hull
     of the rays, Rays.
"""
positive_hull(Rays::Union{Oscar.MatElem,AbstractMatrix}) = Cone(Rays)


"""
    pm_polytope(C::Cone)

Get the underlying polymake `Cone`.
"""
pm_cone(C::Cone) = C.pm_cone

function Base.show(io::IO, C::Cone)
    print(io,"A polyhedral cone")
end
"""
   dim(C::Cone)

Returns the dimension of a cone.
"""
dim(C::Cone) = C.pm_cone.CONE_DIM

"""
   ambient_dim(C::Cone)

Returns the ambient dimension of a cone.
"""
ambient_dim(C::Cone) = C.pm_cone.CONE_AMBIENT_DIM

"""
   codim(C::Cone)

Returns the codimension of a cone.
"""
codim(C::Cone) = ambient_dim(C)-dim(C)

"""
   rays(C::Cone)

Returns the rays of a cone.
"""
rays(C::Cone) = C.pm_cone.RAYS

"""
   facets(C::Cone)

Returns the facets of a cone.
"""
facets(C::Cone) = C.pm_cone.facets


"""
   lineality_space(C::Cone)

   Returns a basis of the lineality space of a cone.
"""
lineality_space(C::Cone) = C.pm_cone.LINEALITY_SPACE

Polymake.visual(C::Cone; opts...) = Polymake.visual(pm_cone(C); opts...)
