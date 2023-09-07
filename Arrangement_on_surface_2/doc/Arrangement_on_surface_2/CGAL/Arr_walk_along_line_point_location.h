
namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2PointLocation

\anchor arr_refwalk_pl

The `Arr_walk_along_line_point_location` class implements a very simple point-location (and
vertical ray-shooting) strategy that improves the naive one.
The algorithm considers an imaginary vertical ray emanating from the
query point, and simulates a walk along the zone of this ray, starting
from the unbounded face until reaching the query point.
In dense arrangements this walk can considerably reduce the number
of traversed arrangement edges, with respect to the na&iuml;ve
algorithm.

The walk-along-a-line point-location object (just like the na&iuml;ve one)
does not use any auxiliary data structures. Thus, attaching it to an
existing arrangement takes constant time, and any ongoing updates to
this arrangement do not affect the point-location object.
It is therefore recommended to use the "walk" point-location strategy
for arrangements that are constantly changing, especially if the number
of issued queries is not large.

\cgalModels{ArrangementPointLocation_2,ArrangementVerticalRayShoot_2}

\sa `ArrangementPointLocation_2`
\sa `ArrangementVerticalRayShoot_2`
\sa `CGAL::Arr_point_location_result<Arrangement>`

*/
template< typename Arrangement >
class Arr_walk_along_line_point_location {
public:

}; /* end Arr_walk_along_line_point_location */
} /* end namespace CGAL */
