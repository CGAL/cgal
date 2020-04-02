
namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3Classes

The class `Surface_mesh_traits_generator_3` provides a type `type`,
that is a model of the concept `SurfaceMeshTraits_3` for the surface
type `Surface`.

\tparam Surface must be a model of the concept `Surface_3`, which means
that it is copy constructible and
assignable. In addition, a `Surface` type is required
<UL>
<LI>either to provide a nested type `Surface::Surface_mesher_traits_3` that is a
model of `SurfaceMeshTraits_3`
<LI>or to be a surface type for which
a specialization of the traits generator
`Surface_mesh_traits_generator_3` exists.
</UL>

Currently, the
library provides partial specializations of the traits generator for
implicit surfaces (`Implicit_surface_3<Traits, Function>`) and gray
level images (`Gray_level_image_3<FT, Point>`).

\sa `SurfaceMeshTraits_3`
\sa `make_surface_mesh`

*/
template< typename Surface >
struct Surface_mesh_traits_generator_3 {

/*!
A model of the concept `SurfaceMeshTraits_3`.
*/
typedef unspecified_type type;

}; /* end Surface_mesh_traits_generator_3 */
} /* end namespace CGAL */
