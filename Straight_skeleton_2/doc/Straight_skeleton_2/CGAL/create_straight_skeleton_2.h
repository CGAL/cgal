namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2

The function `create_exterior_straight_skeleton_2` creates a straight skeleton in the exterior of a 2D polygon with holes. 

The function returns a new `Straight_skeleton_2<K>` in the <I>limited exterior</I> of the 2D polygon \f$ P\f$ given by the point sequence `[vertices_begin,vertices_end]`.
The skeleton in the <I>limited exterior</I> of \f$ P\f$ is the skeleton in the interior of a polygon \f$ Q\f$ with \f$ P\f$ as its hole and a rectangular frame \f$ F\f$ as outer boundary.
The outer boundary \f$ F\f$ is constructed by enlarging the bounding box of \f$ P\f$ a distance \f$ d\f$. 
\f$ d\f$ is a margin sufficiently large to allow an outer offset at dinstance `max_offset` to be obtained from this exterior skeleton, as computed by the function `compute_outer_frame_margin` 


Requirements 
-------------- 

<OL> 
<LI>`K` is any \cgal kernel. 
<LI>`PointIterator::value_type` is equivalent to `K2::Point_2`. 
A cartesian converter is used to convert from `K2::Point_2` to `K::Point_2` 
<LI>`Polygon` is `Polygon_2<K>` or a standard container of `K2::Point_2` elements 
</OL> 

\sa `create_interior_straight_skeleton_2` 
\sa `Straight_skeleton_builder_2<Gt,Ss,Visitor>` 


*/
template<class FT, class PointIterator, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_exterior_straight_skeleton_2 
( FT max_offset
, PointIterator vertices_begin
, PointIterator vertices_end
, K const& k = Exact_predicates_inexact_constructions_kernel
) ;

/*!
\ingroup PkgStraightSkeleton2

The function `create_exterior_straight_skeleton_2` creates a straight skeleton in the exterior of a 2D polygon with holes. 
The function returns a new `Straight_skeleton_2<K>` in the <I>limited exterior</I> of the 2D polygon \f$ P\f$ .
The skeleton in the <I>limited exterior</I> of \f$ P\f$ is the skeleton in the interior of a polygon \f$ Q\f$ with \f$ P\f$ as its hole and a rectangular frame \f$ F\f$ as outer boundary.
The outer boundary \f$ F\f$ is constructed by enlarging the bounding box of \f$ P\f$ a distance \f$ d\f$. 
\f$ d\f$ is a margin sufficiently large to allow an outer offset at dinstance `max_offset` to be obtained from this exterior skeleton, as computed by the function `compute_outer_frame_margin` 
Requirements 
-------------- 

<OL> 
<LI>`K` is any \cgal kernel. 
<LI>`PointIterator::value_type` is equivalent to `K2::Point_2`. 
A cartesian converter is used to convert from `K2::Point_2` to `K::Point_2` 
<LI>`Polygon` is `Polygon_2<K>` or a standard container of `K2::Point_2` elements 
</OL> 

\sa `create_interior_straight_skeleton_2` 
\sa `Straight_skeleton_builder_2<Gt,Ss,Visitor>` 



*/
template<class FT, class Polygon, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_exterior_straight_skeleton_2 ( FT max_offset
, Polygon P
, K const& k = Exact_predicates_inexact_constructions_kernel
) ;

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2

The function `create_offset_polygons_2` creates a straight skeleton in the interior of a 2D polygon with holes. 
The function returns a container with all the offset polygons at distance \f$ offset\f$ obtained from 
the straight skeleton \f$ s\f$. 


<OL> 
<LI>`K` is any \cgal kernel. 
<LI>`FT` is any number type implicitly convertible to `K::FT`. 
<LI>`Straight_skeleton` is `Straight_skeleton_2<K2>`. 
If \f$ K != K2\f$ the straight skeleton is converted to `Straight_skeleton_2<K>`. 
<LI>`Polygon` is a model of `VertexContainer_2`. 
If this first template parameter is omitted, `Polygon_2` is used. 
</OL> 

\sa `create_interior_skeleton_2` 
\sa `create_exterior_skeleton_2` 
\sa `create_interior_skeleton_and_offset_polygons_2` 
\sa `create_exterior_skeleton_and_offset_polygons_2` 


*/
template<class Polygon, class FT, class Straight_skeleton, class K>
std::vector< boost::shared_ptr<Polygon> >
create_offset_polygons_2 ( FT offset
, Straight_skeleton const& s,
, K const& k = Exact_predicates_inexact_constructions_kernel
) ;

/*!
\ingroup PkgStraightSkeleton2

returns a new `Straight_skeleton_2<K>` in the interior of the 2D
polygon with holes whose outer boundary is given by the point sequence
`[outer_contour_vertices_begin,outer_contour_vertices_end]` and its
holes given by `[holes_begin,holes_end]`.

1. `K` is any \cgal kernel.
2. `PointIterator::value_type` is equivalent to `K2::Point_2`.
   A cartesian converter is used to convert from `K2::Point_2` to `K::Point_2`
3. `HoleIterator::value_type` and `Polygon` are `Polygon_2<K>`
   or a standard container of `K2::Point_2` elements 

\sa `create_exterior_straight_skeleton_2`
\sa `Straight_skeleton_builder_2<Gt,Ss>` 

*/
template<class PointIterator, class HoleIterator, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_interior_straight_skeleton_2 ( PointIterator outer_contour_vertices_begin
                                    , PointIterator outer_contour_vertices_end
                                    , HoleIterator  holes_begin
                                    , HoleIterator  holes_end
                                    , K const&      k = Exact_predicates_inexact_constructions_kernel
                                    ) ;

/*!
\ingroup PkgStraightSkeleton2

returns a new `Straight_skeleton_2<K>` in the interior of the 2D
polygon whose outer boundary is given by the point sequence
`[outer_contour_vertices_begin,outer_contour_vertices_end]`.

1. `K` is any \cgal kernel.
2. `PointIterator::value_type` is equivalent to `K2::Point_2`.
   A cartesian converter is used to convert from `K2::Point_2` to `K::Point_2`
3. `HoleIterator::value_type` and `Polygon` are `Polygon_2<K>`
   or a standard container of `K2::Point_2` elements 

\sa `create_exterior_straight_skeleton_2`
\sa `Straight_skeleton_builder_2<Gt,Ss>` 

*/
template<class PointIterator, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_interior_straight_skeleton_2 ( PointIterator outer_contour_vertices_begin
                                    , PointIterator outer_contour_vertices_end
                                    , K const&      k = Exact_predicates_inexact_constructions_kernel
                                    ) ;

/*!
\ingroup PkgStraightSkeleton2

returns a new `Straight_skeleton_2<K>` in the interior of the 2D
polygon `outer_contour`.

1. `K` is any \cgal kernel.
2. `PointIterator::value_type` is equivalent to `K2::Point_2`.
   A cartesian converter is used to convert from `K2::Point_2` to `K::Point_2`
3. `HoleIterator::value_type` and `Polygon` are `Polygon_2<K>`
   or a standard container of `K2::Point_2` elements 

\sa `create_exterior_straight_skeleton_2`
\sa `Straight_skeleton_builder_2<Gt,Ss>` 

*/
template<class Polygon, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_interior_straight_skeleton_2 ( Polygon const& outer_contour
                                    , K const&       k = Exact_predicates_inexact_constructions_kernel
                                    ) ;

} /* namespace CGAL */
