/*!
\ingroup PkgPointSet2Concepts
\cgalConcept

A point set traits class has to provide some primitives that are used by the point set class.
The following catalog lists the involved primitives.
For details about these types see the %Kernel traits documentation.

*/

class PointSetTraits {
public:

/// \name Types
/// @{

/*!

*/
typedef unspecified_type Point_2;

/*!

*/
typedef unspecified_type Circle_2;

/*!

*/
typedef unspecified_type Segment_2;

/*!

*/
typedef unspecified_type FT;

/*!

*/
typedef unspecified_type Orientation_2;

/*!

*/
typedef unspecified_type Side_of_oriented_circle_2;

/*!

*/
typedef unspecified_type Construct_circle_2;

/*!

*/
typedef unspecified_type Compute_squared_distance_2;

/*!

*/
typedef unspecified_type Bounded_side_2;

/*!

*/
typedef unspecified_type Compare_distance_2;

/*!

*/
typedef unspecified_type Construct_center_2;

/// @}

}; /* end PointSetTraits */
