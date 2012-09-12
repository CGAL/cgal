namespace CGAL {

/*!
\ingroup PkgArrangement2

The function `locate` performs a batched point-location operation on a 
given arrangement. It accepts a range of query points, and locates each 
point in the arrangement. The query results are returned through the output 
iterator. Each result is given as a pair of the query point and an object 
representing the arrangement feature that contains it, namely an 
`Object` that may be either `Face_const_handle`, 
`Halfedge_const_handle`, or `Vertex_const_hanlde`. The resulting 
pairs in the output sequence are sorted in increasing \f$ xy\f$-lexicographical 
order of the query points. The function returns a past-the-end iterator of 
the output sequence. 

Requirements 
-------------- 

<UL> 
<LI>`InputIterator::value_type` must be `Traits::Point_2`. 
<LI>`OutputIterator::value_type` must be 
`std::pair<Traits::Point_2,Object>`. 
</UL> 

*/
template<typename Traits, typename Dcel,
typename PointsIterator, typename OutputIterator>
OutputIterator locate (const Arrangement_2<Traits,Dcel>& arr,
PointsIterator points_begin,
PointsIterator points_end,
OutputIterator oi);

} /* namespace CGAL */

