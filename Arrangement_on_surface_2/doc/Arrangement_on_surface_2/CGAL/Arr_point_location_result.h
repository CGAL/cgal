namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2PointLocation
 *
 * A unary metafunction to determine the return type of a point-location or
 * vertical ray-shoot query.
 *
 * \tparam Arrangement must be an instance of the
 * `CGAL::Arrangement_on_surface_2<GeometryTraits,Topology>` class template.
 *
 * \sa `ArrangementPointLocation_2`
 * \sa `ArrangementVerticalRayShoot_2`
 * \sa `CGAL::Arr_naive_point_location<Arrangement>`
 * \sa `CGAL::Arr_walk_along_line_point_location<Arrangement>`
 * \sa `CGAL::Arr_landmarks_point_location<Arrangement,Generator>`
 * \sa `CGAL::Arr_trapezoid_ric_point_location<Arrangement>`
 */
template <typename Arrangement>
struct Arr_point_location_result
{
  /*! The type of the arrangement feature that is the result of a
   * point-location query or a vertical ray-shoot query, namely,
   * `std::variant<Arrangement_on_surface_2::Vertex_const_handle, Arrangement_on_surface_2::Halfedge_const_handle, Arrangement_on_surface_2::Face_const_handle>`
   */
  typedef unspecified_type Type;
}; /* end Arr_point_location_result */

} /* end namespace CGAL */
