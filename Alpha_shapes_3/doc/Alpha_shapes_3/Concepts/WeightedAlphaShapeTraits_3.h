
/*!
\ingroup PkgAlphaShapes3Concepts
\cgalConcept

The concept `WeightedAlphaShapeTraits_3` describes the requirements
for the geometric traits class
of the underlying regular triangulation of a weighted alpha shape.

\cgalRefines `RegularTriangulationTraits_3`, if the underlying triangulation of the alpha shape is a regular triangulation.
\cgalRefines `Periodic_3RegularTriangulationTraits_3`, if the underlying triangulation of the alpha shape is a periodic regular triangulation.

\cgalHasModel All models of `Kernel`.

\sa `CGAL::Exact_predicates_inexact_constructions_kernel` (recommended kernel)
*/

class WeightedAlphaShapeTraits_3 {
public:

/// \name Types
/// @{

/*!
A number type compatible with the type used for
the points coordinates.
*/
typedef unspecified_type FT;

/*!
An object constructor able to compute the squared radius of the
smallest sphere orthogonal to four weighted points `p0, p1, p2, p3`,
and the squared radius of the
smallest sphere orthogonal to three weighted points `p0, p1, p2`,
and the squared radius of smallest sphere orthogonal to
two weighted points `p0, p1`,
and the squared radius of the smallest sphere orthogonal to a single
point `p0`.
*/
typedef unspecified_type Compute_squared_radius_smallest_orthogonal_sphere_3;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
WeightedAlphaShapeTraits_3();

/// @}

/// \name Access Functions
/// @{

/*!

*/
Compute_squared_radius_smallest_orthogonal_sphere_3 compute_squared_radius_smallest_orthogonal_sphere_3_object();

/// @}

}; /* end WeightedAlphaShapeTraits_3 */

