namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Functions

The function `arrange_offset_polygons_2()` arranges the sequence of `Polygon_2` objects obtained by `create_offset_polygons_2()` into `Polygon_with_holes_2` objects by determining geometric parent-hole relationships using a simple algorithm based on the particular characteristics of offset polygons. 

The function determines parent-hole relationships among the polygons given by `[begin,end]` creating 
`boost::shared_ptr< Polygon_with_holes_2<K> >` objects added to the output sequence given `out`.
A `CLOCKWISE` oriented polygon `H` is a hole of a `COUNTERCLOCKWISE` polygon `P`, iff at least one vertex of `H` is `ON_BOUNDED_SIDE` of `P`.

This function should not be used to arrange arbitrary polygons into polygons with holes unless they meet the requirements specified below. 


\tparam K must be a \cgal kernel. 
\tparam InputPolygonPtrIterator must be an input iterator whose `value_type` is a smart pointer 
(such as `boost::shared_ptr`) whose `element_type` is `Polygon_2<K>`. 
The input polygons must be simple. 
The set of input polygons are unique and interior disjoint. That is, given distinct polygons 
`P` and `Q`, there are vertices of `P` which are not on the boundary of `Q` and are all on the 
bounded or unbounded side of `Q` (but not both). 
\tparam OutputPolygonWithHolesPtrIterator must be an output iterator whose `value_type` is a smart pointer 
(such as `boost::shared_ptr`) whose `element_type` is `Polygon_with_holes_2<K>`.  

\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2` 

*/
template<class K, class InputPolygonPtrIterator, class OutputPolygonWithHolesPtrIterator>
void arrange_offset_polygons_2 ( InputPolygonPtrIterator begin
, InputPolygonPtrIterator end
, OutputPolygonWithHolesPtrIterator out
, K const& k 
) ;

} /* namespace CGAL */
