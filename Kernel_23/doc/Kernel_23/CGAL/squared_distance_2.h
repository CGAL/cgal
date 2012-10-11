namespace CGAL {

/*!
\addtogroup squared_distance squared_distance
\ingroup kernel_global_function
\code
#include <squared_distance_2.h> //for 2D functions
#include <squared_distance_3.h> //for 3D functions
\endcode  

*/
/// @{

/*!
computes the square of the Euclidean distance between two geometric
objects.  For arbitrary geometric objects `obj1` and `obj2` the
squared distance is defined as the minimal `squared_distance(p1, p2)`,
where `p1` is a point of `obj1` and `p2` is a point of `obj2`.  Note
that for objects that have an inside (a bounded region), this inside
is part of the object. So, the squared distance from a point inside is
zero, not the squared distance to the closest point on the boundary.

In 2D, the types `Type1` and `Type2` can be any of the following:

- `Point_2`
- `Line_2`
- `Ray_2`
- `Segment_2`
- `Triangle_2`


In 3D, the types `Type1` and `Type2` can be any of the
following:

- `Point_3`
- `Line_3`
- `Ray_3`
- `Segment_3`
- `Plane_3`

\sa `CGAL::compare_distance_to_point` 
\sa `CGAL::compare_signed_distance_to_line` 
\sa `CGAL::compare_signed_distance_to_plane` 
\sa `CGAL::has_larger_distance_to_point` 
\sa `CGAL::has_larger_signed_distance_to_line` 
\sa `CGAL::has_larger_signed_distance_to_plane` 
\sa `CGAL::has_smaller_distance_to_point` 
\sa `CGAL::has_smaller_signed_distance_to_line`
\sa `CGAL::has_smaller_signed_distance_to_plane` 
*/
template <typename Kernel>
Kernel::FT squared_distance(Type1<Kernel> obj1, Type2<Kernel> obj2);
/// @}
}
