
/*!
\ingroup PkgPolyhedronConcepts
\cgalConcept

Required types and member functions for the `PolyhedronTraits_3` concept. This
geometric traits concept is used in the polyhedral surface data
structure `CGAL::Polyhedron_3<Traits>`. Any \cgal kernel is a model of this concept
and can be used directly as template argument.

\cgalRefines `CopyConstructible`
\cgalRefines `Assignable`

\cgalHasModel `CGAL::Polyhedron_traits_3<Kernel>`
\cgalHasModel `CGAL::Polyhedron_traits_with_normals_3<Kernel>`
\cgalHasModel All models of the `Kernel` concept, e.g., `CGAL::Simple_cartesian<FieldNumberType>`

\sa `CGAL::Polyhedron_3<Traits>`

*/

class PolyhedronTraits_3 {
public:

/// \name Types
/// @{

/*!
point type.
*/
typedef unspecified_type Point_3;

/*!
plane equation. Even if plane equations
are not supported with a particular polyhedral surface this
type has to be defined (some dummy type).
*/
typedef unspecified_type Plane_3;

/*!
is an unary function object
that reverses the plane orientation. Must provide `Plane_3 operator()(Plane_3 plane)` that returns the reversed plane. Required
only if plane equations are supported and the `inside_out()`
method is used to reverse the polyhedral surface orientation.
*/
typedef unspecified_type Construct_opposite_plane_3;

/// @}

/// \name Creation
/// @{

/*!

copy constructor.
*/
PolyhedronTraits_3( const PolyhedronTraits_3& traits2);

/*!

assignment.
*/
PolyhedronTraits_3& operator= ( const PolyhedronTraits_3& traits2);

/// @}

/// \name Operations
/// @{

/*!

returns an instance of this function object.
*/
Construct_opposite_plane_3 construct_opposite_plane_3_object();

/// @}

}; /* end PolyhedronTraits_3 */

