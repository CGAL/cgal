
namespace CGAL {

/*!
\ingroup PkgPolyhedronRef

The class `Polyhedron_traits_3` is a model of the `PolyhedronTraits_3`
concept. It defines the geometric types and primitive operations used
in the polyhedral surface data structure
`Polyhedron_3<PolyhedronTraits_3>` in terms of the \cgal `Kernel`. It keeps a local copy of the kernel which makes
it suitable for kernels with local state.

\cgalModels{PolyhedronTraits_3}

\sa `CGAL::Polyhedron_traits_with_normals_3<Kernel>`

\cgalHeading{Implementation}

Since the `PolyhedronTraits_3` concept is a subset of the
concept Kernel, this class just forwards the relevant types and access
member functions from its template argument. However, it is useful
for testing sufficiency of requirements.

\cgalHeading{Example}

Instantiation of a polyhedral surface with the
%Cartesian kernel based on double coordinates.

\cgalExample{Polyhedron/polyhedron_prog_simple.cpp}

*/
template< typename Kernel >
class Polyhedron_traits_3 {
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
typedef Kernel::Plane_3 Plane_3;

/*!

*/
typedef Kernel::Construct_opposite_plane_3
Construct_opposite_plane_3;

/// @}

/// \name Creation
/// @{

/*!
default constructor, uses
`Kernel()` as local reference to the kernel.
*/
Polyhedron_traits_3();

/*!
stores `kernel`
as local reference.
*/
Polyhedron_traits_3( const Kernel& kernel);

/// @}

/// \name Operations
/// @{

/*!

forwarded to `kernel`.
*/
Construct_opposite_plane_3 construct_opposite_plane_3_object();

/// @}

}; /* end Polyhedron_traits_3 */
} /* end namespace CGAL */
