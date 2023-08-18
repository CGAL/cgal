namespace CGAL {

// ---------------------------------------------- INTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2WeightedSkeletonFunctions

\brief creates a weighted straight skeleton in the interior of a 2D polygon with holes.

The outer contour is given by the point sequence `[outer_contour_vertices_begin, outer_contour_vertices_end]`
and its holes are given by `[holes_begin,holes_end[`. Weights of the outer contour are given by
`[outer_contour_weights_begin, outer_contour_weights_end[`, and weights of the holes are given by
`[holes_weights_begin, holes_weights_end]`, in the same order as holes appear in the iterator range.
Within each weight range, weights are given in the same order as the vertices of the contour:
the `i`-th weight in the range is associated to the contour edge between the `i-1`-th and `i`-th vertices.

\tparam PointIterator must be a model of `InputIterator` with value type `InK::Point_2`.
\tparam HoleIterator must be a model of `InputIterator` with `value_type` a model of `SequenceContainer`
                     with value type `InK::Point_2`.
\tparam WeightIterator must be a model of `InputIterator` with `value_type` `InK::FT`.
\tparam HoleWeightsIterator must be a model of `InputIterator` with `value_type` a model of `SequenceContainer`
                            with value type `InK::FT`.
\tparam SsK must be a model of `Kernel`.

\note `Cartesian_converter` and `NT_converter` are used to convert objects from `InK` to `SsK`,
      if they differ.

\pre the range `[outer_contour_vertices_begin, outer_contour_vertices_end[` describes a weakly simple polygon
     that is oriented counterclockwise.
\\pre the range `[holes_begin, holes_end[` describes a sequence of weakly simple polygons that are oriented clockwise.
\pre Holes neither intersect each other nor the outer boundary.
\pre All the weights must be (strictly) positive.
\pre Collinear consecutive contour edges must have equal weights.

\sa `CGAL::create_exterior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2`
*/
template <typename PointIterator, typename HoleIterator,
          typename WeightIterator, typename HoleWeightsIterator,
          typename SsK>
std::shared_ptr< Straight_skeleton_2<SsK> >
create_interior_weighted_straight_skeleton_2(PointIterator outer_contour_vertices_begin,
                                             PointIterator outer_contour_vertices_end,
                                             HoleIterator holes_begin,
                                             HoleIterator holes_end,
                                             WeightIterator outer_contour_weights_begin,
                                             WeightIterator outer_contour_weights_end,
                                             HoleWeightsIterator holes_weights_begin,
                                             HoleWeightsIterator holes_weights_end,
                                             SsK k = CGAL::Exact_predicates_inexact_constructions_kernel());

/*!
\ingroup PkgStraightSkeleton2WeightedSkeletonFunctions

\brief creates a weighted straight skeleton in the interior of a 2D polygon *without* holes.

The outer contour is given by the point sequence `[outer_contour_vertices_begin, outer_contour_vertices_end]`.
Weights of the outer contour are given by `[outer_contour_weights_begin, outer_contour_weights_end]`,
appearing in the same order as the vertices of the contour: the `i`-th weight in the range is associated
to the contour edge between the `i-1`-th and `i`-th vertices.

\tparam PointIterator must be a model of `InputIterator` with value type `InK::Point_2`.
\tparam WeightIterator must be a model of `InputIterator` with `value_type` `InK::FT`.
\tparam SsK must be a model of `Kernel`.

\note `Cartesian_converter` and `NT_converter` are used to convert objects from `InK` to `SsK`,
      if they differ.

\pre the range `[outer_contour_vertices_begin, outer_contour_vertices_end[` describes a weakly simple polygon
     that is oriented counterclockwise.
\pre All the weights must be (strictly) positive.
\pre Collinear consecutive contour edges must have equal weights.

\sa `CGAL::create_exterior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2`
*/
template <typename PointIterator, typename WeightIterator, typename SsK>
std::shared_ptr< Straight_skeleton_2<SsK> >
create_interior_weighted_straight_skeleton_2(PointIterator outer_contour_vertices_begin,
                                             PointIterator outer_contour_vertices_end,
                                             WeightIterator outer_contour_weights_begin,
                                             WeightIterator outer_contour_weights_end,
                                             SsK k = CGAL::Exact_predicates_inexact_constructions_kernel());

/*!
\ingroup PkgStraightSkeleton2WeightedSkeletonFunctions

\brief creates a weighted straight skeleton in the interior of a 2D polygon, possibly with holes.

Range of weights `weights` must be provided in the same order as the contours (i.e., first
the weights of the outer boundary, and then the weights of the holes, if there are any).
Within each range of weights, the weights must be given in the same order as the vertices of the contour:
the `i`-th weight in the range is associated to the contour edge between the `i-1`-th and `i`-th vertices.

\tparam InKPolygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`),
                   or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).
\tparam InKWeights must be a model of `Range` whose value type is itself a model of `Range` with value type `InK::FT`.
\tparam SsK must be a model of `Kernel`.

\note `Cartesian_converter` and `NT_converter` are used to convert objects from `InK` to `SsK`,
      if they differ.

\pre `polygon` is a weakly simple, counterclockwise polygon with clockwise oriented holes.
\pre Holes neither intersect each other nor the outer boundary.
\pre All the weights must be (strictly) positive.
\pre Collinear consecutive contour edges must have equal weights.

\sa `CGAL::create_exterior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2`
*/
template <typename InKPolygon, typename InKWeights, typename SsK>
std::shared_ptr< Straight_skeleton_2<SsK> >
create_interior_weighted_straight_skeleton_2(const InKPolygon& polygon,
                                             const InKWeights& weights,
                                             SsK k = CGAL::Exact_predicates_inexact_constructions_kernel());

// ---------------------------------------------- EXTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2WeightedSkeletonFunctions

\brief creates a weighted straight skeleton in the <I>limited exterior</I> of the 2D polygon `P`
given by the point sequence `[vertices_begin,vertices_end[`.

The skeleton in the <I>limited exterior</I> of `P` is the skeleton in the interior of a polygon `Q`
with `P` as its hole and a rectangular frame `F` as outer boundary.
The outer boundary `F` is constructed by enlarging the bounding box of `P` by a distance `d`.
The distance `d` is a margin sufficiently large to allow an outer offset at distance `max_offset`
to be obtained from this exterior skeleton, as computed by the function `compute_outer_frame_margin()`.

Weights must be provided in the same order as the vertices of the polygon: the `i`-th weight in the range
is associated to the contour edge between the `i-1`-th and `i`-th vertices.

\tparam SsK must be a model of `Kernel`.
\tparam FT must be a model of `FieldNumberType` convertible to `SsK::FT`.
\tparam PointIterator must be a model of `InputIterator` with value type `InK::Point_2`.
\tparam WeightIterator must be a model of `InputIterator` with `value_type` `InK::FT`.

\note `Cartesian_converter` and `NT_converter` are used to convert objects from `InK` to `SsK`,
      if they differ.

\pre `vertices_begin` and `vertices_end` describe a weakly simple polygon that is oriented counterclockwise.
\pre All the weights must be (strictly) positive.
\pre Collinear consecutive contour edges must have equal weights.

\sa `CGAL::create_interior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2`
*/
template <typename FT, typename PointIterator, typename WeightIterator, typename SsK>
std::shared_ptr< Straight_skeleton_2<SsK> >
create_exterior_weighted_straight_skeleton_2(FT max_offset,
                                             PointIterator vertices_begin,
                                             PointIterator vertices_end,
                                             WeightIterator weights_begin,
                                             WeightIterator weights_end,
                                             SsK k = CGAL::Exact_predicates_inexact_constructions_kernel());

/*!
\ingroup PkgStraightSkeleton2WeightedSkeletonFunctions

\brief creates a weighted straight skeleton in the exterior of a 2D polygon (with holes).

The function returns a straight skeleton in the <I>limited exterior</I> of the 2D polygon `P`.
The skeleton in the <I>limited exterior</I> of `P` is the skeleton in the interior of a polygon `Q`
with `P` as its hole and a rectangular frame `F` as outer boundary.
The outer boundary `F` is constructed by enlarging the bounding box of `P` a distance `d`.
`d` is a margin sufficiently large to allow an outer offset at distance `max_offset` to be obtained
from this exterior skeleton, as computed by the function `compute_outer_frame_margin()`.

Weights must be provided in the same order as the vertices of the polygon: the `i`-th weight in the range
is associated to the contour edge between the `i-1`-th and `i`-th vertices.

\warning Holes of the polygon `P` are ignored. If you also need the exterior skeleton for the holes,
         you should call `CGAL::create_interior_straight_skeleton_2()` for each hole.

\tparam SsK must be a model of `Kernel`.
\tparam FT must be a model of `FieldNumberType` convertible to `SsK::FT`.
\tparam InKPolygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                   or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).
\tparam InKWeights must be a model of `Range` whose value type is itself a model of `Range` with value type `InK::FT`.

\note `Cartesian_converter` and `NT_converter` are used to convert objects from `InK` to `SsK`,
      if they differ.

\pre `P` is a weakly simple, counterclockwise polygon with clockwise oriented holes.
\pre Holes neither intersect each other nor the outer boundary.
\pre All the weights must be (strictly) positive.
\pre Collinear consecutive contour edges must have equal weights.
\pre `max_offset` must be positive.

\sa `CGAL::create_interior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2`
*/
template <typename FT, typename Polygon, typename Weights, typename SsK>
std::shared_ptr< Straight_skeleton_2<SsK> >
create_exterior_weighted_straight_skeleton_2(FT max_offset,
                                             const InKPolygon& P,
                                             const InKWeights& weights,
                                             SsK k = CGAL::Exact_predicates_inexact_constructions_kernel());

} /* namespace CGAL */
