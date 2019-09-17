
/*!
\ingroup PkgPolylineDistanceConcepts
\cgalConcept

The concept `PolylineDistanceTraits_2` describes the set of requirements
to be fulfilled to use the function `CGAL::Frechet_distance()`.

\cgalHasModel All models of `Kernel`.
\cgalHasModel `CGAL::Projection_traits_xy_3<K>`
\cgalHasModel `CGAL::Projection_traits_yz_3<K>`
\cgalHasModel `CGAL::Projection_traits_xz_3<K>`

*/

class PolylineDistanceTraits_2 {
public:

/// \name Types
/// @{

/*!
The point type.
*/
typedef unspecified_type Point_2;

/// @}


/// \name Function Objects
/// @{

/*!
A function object to construct a `Segment_2`.

Provides:

`Segment_2 operator()(Point_2 p,Point_2 q)`,

which constructs a segment from two points.
*/
typedef unspecified_type Construct_segment_2;

/*!
A function object to compare the x-coordinate of two points.

Provides the operator:

`bool operator()(Point p, Point q)`

which returns `true` if `p` is before `q`
according to the \f$ x\f$-ordering of points.
*/
typedef unspecified_type Less_x_2;

/*!
A function object to compare the x-coordinate of two points.

Provides the operator:

`Comparison_result operator()(Point p, Point q)`

which returns
`SMALLER, EQUAL` or `LARGER`
according to the
\f$ x\f$-ordering of points `p` and `q`.
*/
typedef unspecified_type Compare_x_2;

/*!
A function object to compare the y-coordinate of two points.
Provides the operator:

`Comparison_result operator()(Point p, Point q)`

which returns
(`SMALLER, EQUAL` or `LARGER`)
according to the
\f$ y\f$-ordering of points `p` and `q`.
*/
typedef unspecified_type Compare_y_2;

/// @}

/// \name Functions
/// @{

/*!

*/
Compare_x_2 compare_x_2_object();

/*!

*/
Compare_y_2 compare_y_2_object();


/// @}

}; /* end PolylineDistanceTraits_2 */

