
namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3Classes

The class `Implicit_surface_3` implements a surface described 
as the zero level 
set of a function \f$ f : \mathbb{R}^3 \longrightarrow \mathbb{R}\f$. 

For this type of surface, the library provides a partial specialization 
of the surface mesher traits generator: 
`Surface_mesh_traits_generator_3<Implicit_surface_3<Traits, Function> >`, 
that provides a traits class, model of the concept 
`SurfaceMeshTraits_3`, 
to be used by the surface mesher. 

\tparam Traits must be a geometric traits class that must be a model of
`ImplicitSurfaceTraits_3`. That concept defines all the types, predicates
and constructors that the `Traits` has to provide to implement the surface
mesh traits
 `Surface_mesh_traits_generator_3<Implicit_surface_3<Traits, Function> >`.


\tparam Function must be a of the concept `ImplicitFunction`. 

The number type `Function::FT` has to match
the type `Traits::FT`. 

\cgalModels `Surface_3`

\sa `make_surface_mesh`
\sa `Surface_mesh_traits_generator_3<Surface>`
\sa `ImplicitSurfaceTraits_3`
\sa `ImplicitFunction`

*/
template< typename Traits, typename Function >
class Implicit_surface_3 {
public:

/// \name Creation 
/// @{

/*!
`f` is the object of type `Function` that represents the implicit 
surface. 

`bounding_sphere` is a bounding sphere of the implicit surface. The 
evaluation of `f` at the center `c` of this sphere must be 
negative: \f$ f(c)<0\f$. 

`error_bound` is a relative error bound 
used to compute intersection points between the implicit surface 
and query segments. This bound is used in the default generated traits class. 
In this traits class, the intersection points between 
the surface and segments/rays/line are constructed by dichotomy. The 
dichotomy is stopped when the size of the intersected 
segment is less than the product of `error_bound` by the 
radius of `bounding_sphere`. 
*/ 
Implicit_surface_3(Function f, 
Sphere_3 bounding_sphere, 
FT error_bound = FT(1e-3)); 

/// @}

}; /* end Implicit_surface_3 */
} /* end namespace CGAL */
