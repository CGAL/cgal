
/*!
\ingroup PkgAlphaShapes2Concepts
\cgalConcept

The concept `WeightedAlphaShapeTraits_2` describes the requirements
for the geometric traits class
of the underlying regular triangulation of a weighted alpha shape.

\cgalRefines `RegularTriangulationTraits_2`, if the underlying triangulation of the alpha shape is a regular triangulation.

\cgalHasModel All models of `Kernel`.
\cgalHasModel Projection traits such as `CGAL::Projection_traits_xy_3<K>`.

\sa `CGAL::Exact_predicates_inexact_constructions_kernel` (recommended kernel)
*/

class WeightedAlphaShapeTraits_2 {
public:

/// \name Types
/// @{

/*!
A coordinate type.
The type must provide a copy constructor, assignment, comparison
operators, negation, multiplication, division and allow the
declaration and initialization with a small integer constant
(cf. requirements for number types).
 An obvious choice would be coordinate type of the point class.
*/
typedef unspecified_type FT;

/// @}

/// \name Creation
/// Only a default constructor is required. Note that further constructors can be provided.
/// @{

/*!
A default constructor.
*/
  WeightedAlphaShapeTraits_2();

/// @}

/// \name Constructions by function objects
/// @{

/*!
Returns an object, which has to be able to compute the squared radius of the
orthogonal circle of the points `p0, p1, p2` or the squared radius of the
smallest orthogonal circle of the points `p0, p1`, as `FT`.
*/
Compute_squared_radius_smallest_orthogonal_circle_2
compute_squared_radius_smallest_orthogonal_circle_2_object();

/// @}

/// \name Predicate by function object
/// @{

/*!
Returns an object, which has to be able to compute the relative position of the
point `test` to the smallest orthogonal circle of the points `p0, p1`.
*/
Power_side_of_bounded_power_circle_2
power_side_of_bounded_power_circle_2_object();

/// @}

}; /* end WeightedAlphaShapeTraits_2 */

