namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Functions

The function `arrange_offset_polygons_2()` arranges the sequence of `Polygon_2` objects obtained by `create_offset_polygons_2()` into `Polygon_with_holes_2` objects by determining geometric parent-hole relationships using a simple algorithm based on the particular characteristics of offset polygons. 

The function determines parent-hole relationships among the polygons given by `[begin,end]` creating 
`boost::shared_ptr< Polygon_with_holes_2<K> >` objects added to the output sequence given `out`.
A `CLOCKWISE` oriented polygon \f$ H\f$ is a hole of a `COUNTERCLOCKWISE` polygon \f$ P\f$, iff at least one vertex of \f$ H\f$ is `ON_BOUNDED_SIDE` of \f$ P\f$.

This function should not be used to arrange arbitrary polygons into polygons with holes unless they meet the requirements specified below. 

### Requirements ###

<OL> 
<LI>`K` is any \cgal kernel. 
<LI>`InputPolygonPtrIterator` is an input iterator whose `value_type` is a `smart ptr` 
(such as `boost::shared_ptr`) whose `element_type` is `Polygon_2<K>`. 
<LI>`InputPolygonPtrIterator` is an output iterator whose `value_type` is a `smart ptr` 
(such as `boost::shared_ptr`) whose `element_type` is `Polygon_with_holes_2<K>`. 
<LI>The input polygons must be simple. 
<LI>The set of input polygons are unique and interior disjoint. That is, given distinct polygons 
\f$ P\f$ and \f$ Q\f$, there are vertices of \f$ P\f$ which are not on the boundary of \f$ Q\f$ and are all on the 
bounded or unbounded side of \f$ Q\f$ (but not both). 
</OL> 

\sa `create_exterior_straight_skeleton_2` 
\sa `Straight_skeleton_builder_2<Gt,Ss>` 

*/
template<class K, class InputPolygonPtrIterator, class OutputPolygonWithHolesPtrIterator>
void arrange_offset_polygons_2 ( InputPolygonPtrIterator begin
, InputPolygonPtrIterator end
, OutputPolygonWithHolesPtrIterator out
, K const& k 
) ;

} /* namespace CGAL */
