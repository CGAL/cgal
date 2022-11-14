namespace CGAL {

/*!
 * \ingroup PkgArrangementOnSurface2PointLocation
 *
 * \anchor arr_reftri_pl
 *
 * The `Arr_triangulation_point_location` class template implements a
 * point-location (and vertical ray-shooting) strategy that is based on
 * triangulation. In particular, the algorithm uses a constrained triangulation,
 * provided by the 2D Triangulations package, as a search strcture. Every time
 * the arrangement is modified the constrained triangulation search-structure is
 * reconstructed from scrach, where the edges of the arrangement are set to be
 * the constrained edges of the triangulation. This strategy is inefficient
 * (especially when the number of modifications applied to the arrangement is
 * high) and provided only for educational purposes.
 *
 * \cgalModels `ArrangementPointLocation_2`
 * \cgalModels `ArrangementVerticalRayShoot_2`
 *
 * \sa `ArrangementPointLocation_2`
 * \sa `ArrangementVerticalRayShoot_2`
 * \sa `CGAL::Arr_point_location_result<Arrangement>`
 * \sa `CGAL_ARR_POINT_LOCATION_VERSION`
 */

template <typename Arrangement_>
class Arr_triangulation_point_location : public Arr_observer<Arrangement_>
{}

}
