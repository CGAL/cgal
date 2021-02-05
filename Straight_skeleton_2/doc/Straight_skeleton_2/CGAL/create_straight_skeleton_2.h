namespace CGAL {

// ---------------------------------------------- EXTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

The function `create_exterior_straight_skeleton_2()` creates a straight skeleton in the exterior
of a 2D polygon with holes.

The function returns a new `Straight_skeleton_2<SsK>` in the <I>limited exterior</I> of the 2D polygon `P`
given by the point sequence `[vertices_begin,vertices_end]`.
The skeleton in the <I>limited exterior</I> of `P` is the skeleton in the interior of a polygon `Q`
with `P` as its hole and a rectangular frame `F` as outer boundary.
The outer boundary `F` is constructed by enlarging the bounding box of `P` a distance `d`.
`d` is a margin sufficiently large to allow an outer offset at distance `max_offset` to be obtained
from this exterior skeleton, as computed by the function `compute_outer_frame_margin()`.

\tparam SsK must be a model of `Kernel`.
\tparam FT must be a model of `FieldNumberType` convertible to `SsK::FT`.
\tparam PointIterator must be a model of `InputIterator` with value type `InK::Point_2`.

\note `Cartesian_converter` is used to convert from `InK::Point_2` to `SsK::Point_2`, if they differ.

\sa `create_interior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2<Traits,Ss,Visitor>`
*/
template<class FT, class PointIterator, class SsK>
boost::shared_ptr< Straight_skeleton_2<SsK> >
create_exterior_straight_skeleton_2( FT max_offset,
                                     PointIterator vertices_begin,
                                     PointIterator vertices_end,
                                     SsK k = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

The function `create_exterior_straight_skeleton_2()` creates a straight skeleton in the exterior
of a 2D polygon with holes.

The function returns a new `Straight_skeleton_2<SsK>` in the <I>limited exterior</I> of the 2D polygon `P`.
The skeleton in the <I>limited exterior</I> of `P` is the skeleton in the interior of a polygon `Q`
with `P` as its hole and a rectangular frame `F` as outer boundary.
The outer boundary `F` is constructed by enlarging the bounding box of `P` a distance `d`.
`d` is a margin sufficiently large to allow an outer offset at distance `max_offset` to be obtained
from this exterior skeleton, as computed by the function `compute_outer_frame_margin()`.

\tparam FT must be a model of `FieldNumberType` convertible to `SsK::FT`.
\tparam SsK must be a model of `Kernel`.
\tparam Polygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).

\note `Cartesian_converter` is used to convert from `InK::Point_2` to `SsK::Point_2`, if they differ.

\sa `create_interior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2<Traits,Ss,Visitor>`
*/
template<class FT, class Polygon, class SsK>
boost::shared_ptr< Straight_skeleton_2<SsK> >
create_exterior_straight_skeleton_2 ( FT max_offset,
                                      const Polygon& P,
                                      SsK k = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

// ---------------------------------------------- INTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

returns a new `Straight_skeleton_2<SsK>` in the interior of the 2D
polygon with holes whose outer boundary is given by the point sequence
`[outer_contour_vertices_begin,outer_contour_vertices_end]` and its
holes given by `[holes_begin,holes_end]`.

\tparam SsK must be a model of `Kernel`.
\tparam PointIterator must be a model of `InputIterator` with value type `InK::Point_2`.
\tparam HoleIterator must be a model of `InputIterator` with `value_type` a model of `SequenceContainer`
                     with value type `InK::Point_2`.

\note `Cartesian_converter` is used to convert from `InK::Point_2` to `SsK::Point_2`, if they differ.

\sa `create_exterior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2<Traits,Ss,Visitor>`
*/
template<class PointIterator, class HoleIterator, class SsK>
boost::shared_ptr< Straight_skeleton_2<SsK> >
create_interior_straight_skeleton_2 ( PointIterator outer_contour_vertices_begin,
                                      PointIterator outer_contour_vertices_end,
                                      HoleIterator holes_begin,
                                      HoleIterator holes_end,
                                      SsK k = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

returns a new `Straight_skeleton_2<SsK>` in the interior of the 2D
polygon whose outer boundary is given by the point sequence
`[outer_contour_vertices_begin,outer_contour_vertices_end]`.

\tparam SsK must be a model of `Kernel`.
\tparam PointIterator must be a model of `InputIterator` with value type `InK::Point_2`.

\note `Cartesian_converter` is used to convert from `InK::Point_2` to `SsK::Point_2`, if they differ.

\sa `create_exterior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2<Traits,Ss,Visitor>`
*/
template<class PointIterator, class SsK>
boost::shared_ptr< Straight_skeleton_2<SsK> >
create_interior_straight_skeleton_2 ( PointIterator outer_contour_vertices_begin,
                                      PointIterator outer_contour_vertices_end,
                                      SsK k = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

returns a new `Straight_skeleton_2<SsK>` in the interior of the 2D
polygon `outer_contour`.

\tparam SsK must be a model of `Kernel`.
\tparam Polygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).

\note `Cartesian_converter` is used to convert from `InK::Point_2` to `SsK::Point_2`,
      if they differ.

\sa `create_exterior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2<Traits,Ss,Visitor>`
*/
template<class Polygon, class SsK>
boost::shared_ptr< Straight_skeleton_2<SsK> >
create_interior_straight_skeleton_2 ( const Polygon& outer_contour,
                                      SsK k = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

} /* namespace CGAL */
