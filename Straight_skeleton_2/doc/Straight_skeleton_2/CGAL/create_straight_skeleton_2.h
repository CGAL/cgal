namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Functions

The function `create_exterior_straight_skeleton_2()` creates a straight skeleton in the exterior of a 2D polygon with holes. 

The function returns a new `Straight_skeleton_2<K>` in the <I>limited exterior</I> of the 2D polygon `P` given by the point sequence `[vertices_begin,vertices_end]`.
The skeleton in the <I>limited exterior</I> of `P` is the skeleton in the interior of a polygon `Q` with `P` as its hole and a rectangular frame `F` as outer boundary.
The outer boundary `F` is constructed by enlarging the bounding box of `P` a distance `d`. 
`d` is a margin sufficiently large to allow an outer offset at distance `max_offset` to be obtained from this exterior skeleton, as computed by the function `compute_outer_frame_margin()` 


\cgalHeading{Requirements}

\tparam FT a number type
\tparam K is any \cgal Kernel
\tparam PointIterator an `InputIterator` with `value_type` being equivalent to `K2::Point_2`
`Cartesian_converter` is used to convert from `K2::Point_2` to `K::Point_2`

\sa `create_interior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`


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
\ingroup PkgStraightSkeleton2Functions

The function `create_exterior_straight_skeleton_2()` creates a straight skeleton in the exterior of a 2D polygon with holes. 
The function returns a new `Straight_skeleton_2<K>` in the <I>limited exterior</I> of the 2D polygon `P`.
The skeleton in the <I>limited exterior</I> of `P` is the skeleton in the interior of a polygon `Q` with `P` as its hole and a rectangular frame `F` as outer boundary.
The outer boundary `F` is constructed by enlarging the bounding box of `P` a distance `d`. 
`d` is a margin sufficiently large to allow an outer offset at distance `max_offset` to be obtained from this exterior skeleton, as computed by the function `compute_outer_frame_margin()` 

\cgalHeading{Requirements}

\tparam FT a number type
\tparam K is any \cgal kernel. 
\tparam Polygon is `Polygon_2<K>` or a standard container of `K::Point_2` elements.

\sa `create_interior_straight_skeleton_2()` 
\sa `Straight_skeleton_builder_2` 
*/
template<class FT, class Polygon, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_exterior_straight_skeleton_2 ( FT max_offset
, Polygon P
, K const& k = Exact_predicates_inexact_constructions_kernel
) ;

/*!
\ingroup PkgStraightSkeleton2Functions

returns a new `Straight_skeleton_2<K>` in the interior of the 2D
polygon with holes whose outer boundary is given by the point sequence
`[outer_contour_vertices_begin,outer_contour_vertices_end]` and its
holes given by `[holes_begin,holes_end]`.

\tparam K is any \cgal Kernel
\tparam PointIterator an `InputIterator` with `value_type` begin equivalent to `K2::Point_2`.
        `Cartesian_converter` is used to convert from `K2::Point_2` to `K::Point_2`
\tparam HoleIterator an `InputIterator` with `value_type` being `Polygon_2<K>`
        or a standard container of `K2::Point_2` elements.


\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`

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
\ingroup PkgStraightSkeleton2Functions

returns a new `Straight_skeleton_2<K>` in the interior of the 2D
polygon whose outer boundary is given by the point sequence
`[outer_contour_vertices_begin,outer_contour_vertices_end]`.

\tparam K is any \cgal Kernel
\tparam PointIterator an `InputIterator` with `value_type` begin equivalent to `K2::Point_2`.
        `Cartesian_converter` is used to convert from `K2::Point_2` to `K::Point_2`

\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2` 

*/
template<class PointIterator, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_interior_straight_skeleton_2 ( PointIterator outer_contour_vertices_begin
                                    , PointIterator outer_contour_vertices_end
                                    , K const&      k = Exact_predicates_inexact_constructions_kernel
                                    ) ;

/*!
\ingroup PkgStraightSkeleton2Functions

returns a new `Straight_skeleton_2<K>` in the interior of the 2D
polygon `outer_contour`.

\tparam K is any \cgal Kernel
\tparam Polygon an `InputIterator` with `value_type` being `Polygon_2<K>`
        or a standard container of `K2::Point_2` elements.

\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`

*/
template<class Polygon, class K>
boost::shared_ptr< Straight_skeleton_2<K> >
create_interior_straight_skeleton_2 ( Polygon const& outer_contour
                                    , K const&       k = Exact_predicates_inexact_constructions_kernel
                                    ) ;

} /* namespace CGAL */
