/*!
\ingroup PkgArrangementOnSurface2PointLocation

The macro `CGAL_ARR_POINT_LOCATION_VERSION` can be used to configure
the point-location query API. In particular, it determines which version
of the result type of the point-location and vertical ray-shooting queries
should be used by models of the concepts `ArrangementPointLocation_2`
and `ArrangementVerticalRayShoot_2`, and by the free function
`locate`. The `CGAL_ARR_POINT_LOCATION_VERSION` should be defined before any \cgal header
is included.

- `CGAL_ARR_POINT_LOCATION_VERSION` == 1, the result type is set to be `CGAL::Object`.
- `CGAL_ARR_POINT_LOCATION_VERSION` == 2, the result type is set to be
`boost::variant<Vertex_const_handle,Halfedge_const_handle,Face_const_handle>`, where `Vertex_const_handle`, `Halfedge_const_handle`, and
`Face_const_handle` are the corresponding nested types in a `CGAL::Arrangement_2` instance.


\sa `ArrangementPointLocation_2`
\sa `ArrangementVerticalRayShoot_2`
\sa `CGAL::Arr_point_location_result<Arrangement>`
*/
#define CGAL_ARR_POINT_LOCATION_VERSION

namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2PointLocation

A unary metafunction to determine the return type of a point-location
or vertical ray-shoot query.

\tparam Arrangement must be an instance of the `CGAL::Arrangement_2<Traits,Dcel>` class template.

\sa `ArrangementPointLocation_2`
\sa `ArrangementVerticalRayShoot_2`
\sa `CGAL::Arr_naive_point_location<Arrangement>`
\sa `CGAL::Arr_walk_along_line_point_location<Arrangement>`
\sa `CGAL::Arr_landmarks_point_location<Arrangement,Generator>`
\sa `CGAL::Arr_trapezoid_ric_point_location<Arrangement>`
\sa `CGAL_ARR_POINT_LOCATION_VERSION`
*/
template <class Arrangement>
struct Arr_point_location_result
{
  /*! The type of the arrangement feature that is the result of a
   * point-location query or a vertical ray-shoot query, namely,
   * `boost::variant<Arrangement::Vertex_const_handle, Arrangement::Halfedge_const_handle, Arrangement::Face_const_handle>`
   * if `::CGAL_ARR_POINT_LOCATION_VERSION` == 2, which is the default, otherwise
   * `CGAL::Object`.
   */
  typedef unspecified_type Type;
}; /* end Arr_point_location_result */
} /* end namespace CGAL */
