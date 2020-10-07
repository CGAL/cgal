namespace CGAL {

// ---------------------------------------------- EXTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

The function `create_exterior_straight_skeleton_2()` creates a straight skeleton in the exterior
of a 2D polygon with holes.

The function returns a new `Straight_skeleton_2<SkeletonK>` in the <I>limited exterior</I> of the 2D polygon `P`
given by the point sequence `[vertices_begin,vertices_end]`.
The skeleton in the <I>limited exterior</I> of `P` is the skeleton in the interior of a polygon `Q`
with `P` as its hole and a rectangular frame `F` as outer boundary.
The outer boundary `F` is constructed by enlarging the bounding box of `P` a distance `d`.
`d` is a margin sufficiently large to allow an outer offset at distance `max_offset` to be obtained
from this exterior skeleton, as computed by the function `compute_outer_frame_margin()`

\tparam SkeletonK must be a model of `Kernel`.
\tparam FT must be a model of `FieldType` convertible to `SkeletonK::FT`.
\tparam PointIterator must be a model of `InputIterator` with `value_type` being equivalent to `InputKernel::Point_2`.


\note `Cartesian_converter` is used to convert from `InputKernel::Point_2` to `SkeletonK::Point_2`.

\sa `create_interior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`
*/
template<class FT, class PointIterator, class SkeletonK>
boost::shared_ptr< Straight_skeleton_2<SkeletonK> >
create_exterior_straight_skeleton_2( FT max_offset,
                                     PointIterator vertices_begin,
                                     PointIterator vertices_end,
                                     SkeletonK k = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

The function `create_exterior_straight_skeleton_2()` creates a straight skeleton in the exterior
of a 2D polygon with holes.

The function returns a new `Straight_skeleton_2<SkeletonK>` in the <I>limited exterior</I> of the 2D polygon `P`.
The skeleton in the <I>limited exterior</I> of `P` is the skeleton in the interior of a polygon `Q`
with `P` as its hole and a rectangular frame `F` as outer boundary.
The outer boundary `F` is constructed by enlarging the bounding box of `P` a distance `d`.
`d` is a margin sufficiently large to allow an outer offset at distance `max_offset` to be obtained
from this exterior skeleton, as computed by the function `compute_outer_frame_margin()`

\tparam FT must be a model of `FieldType` convertible to `SkeletonK::FT`.
\tparam SkeletonK is any \cgal kernel.
\tparam Polygon must be a model of `SequenceContainer` with value type `InputKernel::Point_2` elements
                (e.g. `Polygon_2<InputKernel>` or `Polygon_with_holes_2<InputKernel>`).

\note `Cartesian_converter` is used to convert from `InputKernel::Point_2` to `SkeletonK::Point_2`.

\sa `create_interior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`
*/
template<class FT, class Polygon, class SkeletonK>
boost::shared_ptr< Straight_skeleton_2<SkeletonK> >
create_exterior_straight_skeleton_2 ( FT max_offset,
                                      const Polygon& P,
                                      SkeletonK k = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

// ---------------------------------------------- INTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

returns a new `Straight_skeleton_2<SkeletonK>` in the interior of the 2D
polygon with holes whose outer boundary is given by the point sequence
`[outer_contour_vertices_begin,outer_contour_vertices_end]` and its
holes given by `[holes_begin,holes_end]`.

\tparam SkeletonK must be a model of `Kernel`.
\tparam PointIterator must be a model of `InputIterator` with `value_type` being equivalent to `InputKernel::Point_2`.
\tparam HoleIterator must be a model of `InputIterator` with `value_type` a model of `SequenceContainer`
                     with value type `InputKernel::Point_2`.

\note `Cartesian_converter` is used to convert from `InputKernel::Point_2` to `SkeletonK::Point_2`.

\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`
*/
template<class PointIterator, class HoleIterator, class SkeletonK>
boost::shared_ptr< Straight_skeleton_2<SkeletonK> >
create_interior_straight_skeleton_2 ( PointIterator outer_contour_vertices_begin,
                                      PointIterator outer_contour_vertices_end,
                                      HoleIterator holes_begin,
                                      HoleIterator holes_end,
                                      SkeletonK k = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

returns a new `Straight_skeleton_2<SkeletonK>` in the interior of the 2D
polygon whose outer boundary is given by the point sequence
`[outer_contour_vertices_begin,outer_contour_vertices_end]`.

\tparam SkeletonK must be a model of `Kernel`.
\tparam PointIterator must be a model of `InputIterator` with value type `InputKernel::Point_2`.

\note `Cartesian_converter` is used to convert from `InputKernel::Point_2` to `SkeletonK::Point_2`.

\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`
*/
template<class PointIterator, class SkeletonK>
boost::shared_ptr< Straight_skeleton_2<SkeletonK> >
create_interior_straight_skeleton_2 ( PointIterator outer_contour_vertices_begin,
                                      PointIterator outer_contour_vertices_end,
                                      SkeletonK k = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

returns a new `Straight_skeleton_2<SkeletonK>` in the interior of the 2D
polygon `outer_contour`.

\tparam SkeletonK must be a model of `Kernel`.
\tparam Polygon must be a model of `SequenceContainer` with value type `InputKernel::Point_2` elements
                (e.g. `Polygon_2<InputKernel>` or `Polygon_with_holes_2<InputKernel>`).

\note `Cartesian_converter` is used to convert from `InputKernel::Point_2` to `SkeletonK::Point_2`.

\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`
*/
template<class Polygon, class SkeletonK>
boost::shared_ptr< Straight_skeleton_2<SkeletonK> >
create_interior_straight_skeleton_2 ( const Polygon& outer_contour,
                                      SkeletonK k = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

} /* namespace CGAL */
