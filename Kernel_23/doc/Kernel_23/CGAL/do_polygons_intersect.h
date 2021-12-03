namespace CGAL {


/*!
\addtogroup do_intersect_polygons_grp
\ingroup do_intersect

\sa `do_intersect_circular_grp`
\sa `do_intersect_spherical_grp`
\sa `do_intersect_linear_grp`
\sa `intersection_grp`
*/

/// @{
/*!
checks whether the polygons `p1` and `p2` intersect.
`Polygon1` and `Polygon2` can be an instantiation of one of the two
polygon classes, namely `CGAL::Polygon_2<K,C>` or `CGAL::Polygon_with_holes_2<K,C>`.
Two polygons are considered as intersecting iff there exists a point that is in both regions
enclosed by the polygons (boundary included). If you are looking for a regularized
do-intersect function (i.e. existence of an infinitesimal disk covered by both polygons),
see the function \link boolean_do_intersect `CGAL::Regularized_boolean_operations_2::do_intersect()` \endlink.

\note This function depends on the package \ref PkgBooleanSetOperations2 and is documented here for convenience.

*/
template <class Polygon1, class Polygon2>
bool do_intersect(const Polygon1& p1, const Polygon2& p2);

}
