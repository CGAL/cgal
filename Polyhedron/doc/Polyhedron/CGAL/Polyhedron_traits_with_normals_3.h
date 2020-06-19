
namespace CGAL {

/*!
\ingroup PkgPolyhedronRef

The class `Polyhedron_traits_with_normals_3` is a model of the `PolyhedronTraits_3`
concept. It defines the geometric types and primitive operations used
in the polyhedral surface data structure
`Polyhedron_3<PolyhedronTraits_3>`. `Polyhedron_traits_with_normals_3` uses the
normal vector from `Kernel` for the plane equation in facets. It
keeps a local copy of the kernel which makes it suitable for kernels
with local state.

\cgalModels `PolyhedronTraits_3`

\sa `CGAL::Polyhedron_traits_3<Kernel>`

\cgalHeading{Example}

We use this traits class to instantiate a polyhedral surface with a
normal vector and no plane equation for each facet. We compute the
normal vector assuming exact arithmetic (integers in this example)
and convex planar facets.

\cgalExample{Polyhedron/polyhedron_prog_normals.cpp}

*/
template< typename Kernel >
class Polyhedron_traits_with_normals_3 {
public:

/// \name Types
/// @{

/*!
the `Kernel` model.
*/
typedef unspecified_type Kernel;

/*!

*/
typedef Kernel::Point_3 Point_3;

/*!

*/
typedef Kernel::Vector_3 Plane_3;

/*!

*/
typedef Kernel::Construct_opposite_vector_3
Construct_opposite_plane_3;

/// @}

/// \name Creation
/// @{

/*!
default constructor, uses
`Kernel()` as local reference to the kernel.
*/
Polyhedron_traits_with_normals_3();

/*!
stores `kernel` as local reference.
*/
Polyhedron_traits_with_normals_3( const Kernel& kernel);

/// @}

/// \name Operations
/// @{

/*!

forwarded to `kernel`.
*/
Construct_opposite_plane_3 construct_opposite_plane_3_object();

/// @}

}; /* end Polyhedron_traits_with_normals_3 */
} /* end namespace CGAL */
