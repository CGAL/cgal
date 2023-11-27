namespace CGAL {

/*! \ingroup PkgArrangementOnSurface2PointLocation
 *
 * performs a batched point-location operation on a given arrangement. It
 * accepts a collection of query points, locates each point in a given
 * arrangement, and inserts the query results into an output container given
 * through an output iterator. Each query result is given as a pair of the query
 * point and an object representing the arrangement feature that contains it,
 * namely a discriminated union container of the types `Face_const_handle`,
 * `Halfedge_const_handle`, and `Vertex_const_hanlde`. The resulting pairs in
 * the output container are sorted in increasing \f$xy\f$-lexicographical order
 * of the query points.
 *
 * \param arr The arrangement.
 * \param begin The begin iterator of the container of input points.
 * \param end The past-the-end iterator of the container of input points.
 * \param oi The output iterator that points at the output container.
 * \return The past-the-end iterator of the output container.
 *
 * \cgalHeading{Requirements}
 *
 * \pre The value type of `InputIterator` must be convertible to
 * `Arrangement_2::Point_2`.
 * \pre Dereferencing `oi` must yield an  object convertible to
 * `std::pair<Arrangement_2::Point_2, Arr_point_location_result<Arrangement_2>::%Type>`.
 *
 * \sa `CGAL::Arr_point_location_result<Arrangement>`
 *
 */
template <typename Traits, typename Dcel,
typename InputIterator, typename OutputIterator>
OutputIterator locate (const Arrangement_2<Traits, Dcel>& arr,
                       InputIterator begin,
                       InputIterator end,
                       OutputIterator oi);

} /* namespace CGAL */
