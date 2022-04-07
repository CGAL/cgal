
/*!
\ingroup PkgAlphaShapes2Concepts
\cgalConcept

The concept `AlphaShapeTraits_2` describes the requirements for the geometric traits
class of the underlying Delaunay triangulation of a basic alpha shape.

\cgalRefines `DelaunayTriangulationTraits_2`, if the underlying triangulation of the alpha shape is a Delaunay triangulation.
\cgalRefines `Periodic_2DelaunayTriangulationTraits_2`, if the underlying triangulation of the alpha shape is a periodic Delaunay triangulation.

\cgalHasModel All models of `Kernel`.
\cgalHasModel Projection traits such as `CGAL::Projection_traits_xy_3<K>`.

\sa `CGAL::Exact_predicates_inexact_constructions_kernel` (recommended kernel)
*/

class AlphaShapeTraits_2 {
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
  AlphaShapeTraits_2();

/// @}

/// \name Constructions by function objects
/// @{

/*!
Returns an object, which has to be able to compute the squared radius of the
circle of the points `p0, p1, p2` or the squared radius of smallest circle
of the points `p0, p1`, as `FT` associated with <I>the metric used
by `Dt`</I>.
*/
Compute_squared_radius_2 compute_squared_radius_2_object();

/// @}

/// \name Predicate by function object
/// @{

/*!
Returns an object, which has to be able to compute the relative position of
point `test` to the smallest circle of the points `p0, p1`, using
<I>the same metric as `Dt`</I>.
*/
Side_of_bounded_circle_2 side_of_bounded_circle_2_object();

/// @}

}; /* end AlphaShapeTraits_2 */

