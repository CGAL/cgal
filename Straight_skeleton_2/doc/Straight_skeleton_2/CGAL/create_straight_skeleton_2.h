namespace CGAL {

// ---------------------------------------------- INTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

\brief creates a straight skeleton in the interior of the 2D polygon with holes.

The outer boundary of the polygon is given by the point sequence `[outer_contour_vertices_begin,outer_contour_vertices_end[`
and its holes given by `[holes_begin, holes_end[`.

\tparam PointIterator must be a model of `InputIterator` with value type `InK::Point_2`.
\tparam HoleIterator must be a model of `InputIterator` with `value_type` a model of `SequenceContainer`
                     with value type `InK::Point_2`.
\tparam SsK must be a model of `Kernel`.

\note `Cartesian_converter` is used to convert from `InK::Point_2` to `SsK::Point_2`, if they differ.

\pre the range `[outer_contour_vertices_begin, outer_contour_vertices_end[` describes a weakly simple polygon
      that is oriented counterclockwise.
\\pre the range `[holes_begin, holes_end[` describes a sequence of weakly simple polygons that are oriented clockwise.
\pre Holes neither intersect each other nor the outer boundary.

\sa `CGAL::create_exterior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2`
*/
template <typename PointIterator, typename HoleIterator, typename SsK>
std::shared_ptr< Straight_skeleton_2<SsK> >
create_interior_straight_skeleton_2(PointIterator outer_contour_vertices_begin,
                                    PointIterator outer_contour_vertices_end,
                                    HoleIterator holes_begin,
                                    HoleIterator holes_end,
                                    SsK k = CGAL::Exact_predicates_inexact_constructions_kernel());

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

\brief creates a straight skeleton in the interior of the 2D polygon.

The outer boundary of the polygon is given by the point sequence `[outer_contour_vertices_begin,outer_contour_vertices_end[`.

\tparam PointIterator must be a model of `InputIterator` with value type `InK::Point_2`.
\tparam SsK must be a model of `Kernel`.

\note `Cartesian_converter` is used to convert from `InK::Point_2` to `SsK::Point_2`, if they differ.

\pre the range `[outer_contour_vertices_begin, outer_contour_vertices_end[` describes a weakly simple polygon
      that is oriented counterclockwise.

\sa `CGAL::create_exterior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2`
*/
template <typename PointIterator, typename SsK>
std::shared_ptr< Straight_skeleton_2<SsK> >
create_interior_straight_skeleton_2(PointIterator outer_contour_vertices_begin,
                                    PointIterator outer_contour_vertices_end,
                                    SsK k = CGAL::Exact_predicates_inexact_constructions_kernel());

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

\brief creates a straight skeleton in the interior of the 2D polygon `polygon`.

\warning Holes of the polygon are ignored. If you also need the exterior skeleton for the holes,
         you should call `CGAL::create_interior_straight_skeleton_2()` for each hole.

\tparam Polygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).
\tparam SsK must be a model of `Kernel`.

\note `Cartesian_converter` is used to convert from `InK::Point_2` to `SsK::Point_2`,
      if they differ.

\pre `polygon` is a weakly simple, counterclockwise polygon with clockwise oriented holes.
\pre Holes neither intersect each other nor the outer boundary.

\sa `CGAL::create_exterior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2`
*/
template <typename Polygon, typename SsK>
std::shared_ptr< Straight_skeleton_2<SsK> >
create_interior_straight_skeleton_2(const Polygon& polygon,
                                    SsK k = CGAL::Exact_predicates_inexact_constructions_kernel());

// ---------------------------------------------- EXTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

\brief creates a straight skeleton in the exterior of a 2D polygon with holes.

The function returns a straight skeleton in the <I>limited exterior</I> of the 2D polygon `P`
given by the point sequence `[vertices_begin,vertices_end[`.
The skeleton in the <I>limited exterior</I> of `P` is the skeleton in the interior of a polygon `Q`
with `P` as its hole and a rectangular frame `F` as outer boundary.
The outer boundary `F` is constructed by enlarging the bounding box of `P` a distance `d`.
`d` is a margin sufficiently large to allow an outer offset at distance `max_offset` to be obtained
from this exterior skeleton, as computed by the function `compute_outer_frame_margin()`.

\tparam SsK must be a model of `Kernel`.
\tparam FT must be a model of `FieldNumberType` convertible to `SsK::FT`.
\tparam PointIterator must be a model of `InputIterator` with value type `InK::Point_2`.

\note `Cartesian_converter` is used to convert from `InK::Point_2` to `SsK::Point_2`, if they differ.

\pre `vertices_begin` and `vertices_end` describe a weakly simple polygon that is oriented counterclockwise.

\sa `CGAL::create_interior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2`
*/
template <typename FT, typename PointIterator, typename SsK>
std::shared_ptr< Straight_skeleton_2<SsK> >
create_exterior_straight_skeleton_2(FT max_offset,
                                    PointIterator vertices_begin,
                                    PointIterator vertices_end,
                                    SsK k = CGAL::Exact_predicates_inexact_constructions_kernel());

/*!
\ingroup PkgStraightSkeleton2SkeletonFunctions

\brief creates a straight skeleton in the exterior of a 2D polygon with holes.

The function returns a straight skeleton in the <I>limited exterior</I> of the 2D polygon `P`.
The skeleton in the <I>limited exterior</I> of `P` is the skeleton in the interior of a polygon `Q`
with `P` as its hole and a rectangular frame `F` as outer boundary.
The outer boundary `F` is constructed by enlarging the bounding box of `P` a distance `d`.
`d` is a margin sufficiently large to allow an outer offset at distance `max_offset` to be obtained
from this exterior skeleton, as computed by the function `compute_outer_frame_margin()`.

\tparam SsK must be a model of `Kernel`.
\tparam FT must be a model of `FieldNumberType` convertible to `SsK::FT`.
\tparam Polygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).

\note `Cartesian_converter` is used to convert from `InK::Point_2` to `SsK::Point_2`, if they differ.

\pre `P` is a weakly simple, counterclockwise polygon with clockwise oriented holes.
\pre Holes neither intersect each other nor the outer boundary.
\pre All the weights must be (strictly) positive.
\pre `max_offset` is positive.

\sa `CGAL::create_interior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2`
*/
template <typename FT, typename Polygon, typename SsK>
std::shared_ptr< Straight_skeleton_2<SsK> >
create_exterior_straight_skeleton_2(FT max_offset,
                                    const Polygon& P,
                                    SsK k = CGAL::Exact_predicates_inexact_constructions_kernel());

} /* namespace CGAL */
