namespace CGAL {

/*!
\ingroup PkgArrangementOnSurface2PointLocation

Performs a batched point-location operation on a
given arrangement. It accepts a range of query points, and locates each
point in the arrangement. The query results are returned through the output
iterator. Each query result is given as a pair of the query point and an
object representing the arrangement feature that contains it, namely a
discriminated union container of the bounded types `Face_const_handle`,
`Halfedge_const_handle`, and `Vertex_const_hanlde`. The resulting
pairs in the output sequence are sorted in increasing \f$ xy\f$-lexicographical
order of the query points. The function returns a past-the-end iterator of
the output sequence.

\cgalHeading{Requirements}

<UL>
<LI>`InputIterator::value_type` must be `Arrangement_2::Point_2`.
<LI>`OutputIterator::value_type` must be convertible to
`std::pair<Arrangement_2::Point_2, Arr_point_location_result<Arrangement_2>::%Type>`.
</UL>

\cgalHeading{A Note on Backwards Compatibility}
This function used to return `CGAL::Object` up to
\cgal version 4.2. Starting with \cgal version 4.3 the return type
is determined by the metafunction `CGAL::Arr_point_location_result`.
To preserve backwards compatibility
`CGAL::Object` can be constructed from the new return type
implicitly, but switching to the new style is recommended. To enable
the old style without any overhead, the macro
`::CGAL_ARR_POINT_LOCATION_VERSION` can be defined to 1 before any
\cgal header is included.

\sa `CGAL::Arr_point_location_result<Arrangement>`
\sa `CGAL_ARR_POINT_LOCATION_VERSION`

*/
template<typename Traits, typename Dcel,
typename InputIterator, typename OutputIterator>
OutputIterator locate (const Arrangement_2<Traits,Dcel>& arr,
InputIterator points_begin,
InputIterator points_end,
OutputIterator oi);

} /* namespace CGAL */

