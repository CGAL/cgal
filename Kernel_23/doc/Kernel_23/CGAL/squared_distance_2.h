namespace CGAL {

/*!
\addtogroup squared_distance_grp

\code
#include <CGAL/squared_distance_2.h> //for 2D functions
#include <CGAL/squared_distance_3.h> //for 3D functions
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

as well as any combination of `Kernel::Point_2` and `Kernel::Weighted_point_2`

In 3D, the types `Type1` and `Type2` can be any of the
following:

- `Point_3`
- `Line_3`
- `Ray_3`
- `Segment_3`
- `Plane_3`

as well as combinations `Point_3`/`Triangle_3`,
and `Weighted_point_3`/`Triangle_3`.

\sa `compare_distance_to_point_grp`
\sa `compare_signed_distance_to_line_grp`
\sa `compare_signed_distance_to_plane_grp`
\sa `has_larger_distance_to_point_grp`
\sa `has_larger_signed_distance_to_line_grp`
\sa `has_larger_signed_distance_to_plane_grp`
\sa `has_smaller_distance_to_point_grp`
\sa `has_smaller_signed_distance_to_line_grp`
\sa `has_smaller_signed_distance_to_plane_grp`
*/
template <typename Kernel>
Kernel::FT squared_distance(Type1<Kernel> obj1, Type2<Kernel> obj2);

/*!
computes the squared distance and information on where the closest point of `seg` is located.

\code{cpp}
template <typename FT>
struct Squared_distance_dimension_index {
  FT squared_distance;
  unsigned int dimension;
  unsigned int index;
  operator FT() const { return squared_distance;}
};
\endcode

@param pt the input point
@param seg the input segment
@returns an object `res` where `res.squared_distance` is the squared distance, where `res.dimension` is `0`, or `1` if the closest point is on a vertex or on the interior of the segment `seg`, respectively,
and where `res.index` is `0`, or `1`, if `res.dimension==0` and the closest point is on the source or target of the segment, respectively.
 */
template <typename Kernel>
Squared_distance_dimension_index<Kernel::FT> squared_distance(Point_3<Kernel> pt, Segment_3<Kernel> seg, Tag_true);

/*!
computes the squared distance between `pt` and `tri`, and information on where the point of `tri` closest to `pt` is located.

\code{cpp}
template <typename FT>
struct Squared_distance_dimension_index {
  FT squared_distance;
  unsigned int dimension;
  unsigned int index;
  operator FT() const { return squared_distance;}
};
\endcode

@param pt the input point
@param tri the input segment

@returns an object `res` where `res.squared_distance` is the squared distance, where `res.dimension` is `0`, `1`, or`2`, if the closest point is on a vertex, on an edge, or on the interior of the triangle `tri`, respectively,
and where `res.index` is `0`, `1`, or `2`, if `res.dimension==1` and the closest point is on the edge opposite to the ith vertex of the triangle.


 */
template <typename Kernel>
Squared_distance_dimension_index<Kernel::FT> squared_distance(Point_3<Kernel> pt, Triangle_3<Kernel> tri, Tag_true);
/// @}
}
