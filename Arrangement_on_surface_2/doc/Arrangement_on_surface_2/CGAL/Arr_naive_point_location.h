namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2PointLocation
 *
 * \anchor arr_refnaive_pl
 *
 * The `Arr_naive_point_location` class implements a naive algorithm that
 * traverses all the vertices and halfedges in the arrangement in search for an
 * answer to a point-location query.  The query time is therefore linear in the
 * complexity of the arrangement.  Naturally, this point-location strategy could
 * turn into a heavy time-consuming process when applied to dense arrangements.
 *
 * \cgalModels{AosPointLocation_2,AosVerticalRayShoot_2}
 *
 * \sa `AosPointLocation_2`
 * \sa `AosVerticalRayShoot_2`
 * \sa `CGAL::Arr_point_location_result<Arrangement>`
 */
template <typename Arrangement>
class Arr_naive_point_location {};

} /* end namespace CGAL */
