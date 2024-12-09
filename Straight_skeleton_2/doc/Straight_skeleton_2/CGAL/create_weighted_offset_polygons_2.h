namespace CGAL {

// ---------------------------------------------- INTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2WeightedOffsetFunctions

\brief returns a container with all the inner offset polygons at distance `offset`
      of a 2D polygon with holes.

The outer boundary is `outer_boundary` and its holes are given by `[holes_begin,holes_end)`.
Weights of the outer contour are given by `outer_boundary_weights` and weights of the holes
are given by `[holes_weights_begin, holes_weights_end]`, in the same order as holes appear
in the iterator range. Within each weight range, weights are given in the same order as the vertices
of the contour: the `i`-th weight in the range is associated to the contour edge between the `i-1`-th
and `i`-th vertices.

\note A temporary straight skeleton is constructed in the interior of the input polygon to obtain the offsets.
The construction of this skeleton is the most expensive operation, therefore, to construct offsets
at more than one single distance, it is advised to use the separate functions `create_interior_straight_skeleton_2()`
and `create_offset_polygons_2()` instead.

\tparam OfK must be a model of `Kernel`. It is used to instantiate
            `Polygon_offset_builder_traits_2<OfK>` for constructing the offset polygons.
\tparam SsK must be a model of `Kernel`. It is used to instantiate
            `Straight_skeleton_builder_traits_2<SsK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldNumberType` convertible to `OfK::FT` and `SsK::FT`.
\tparam InKPolygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`).
\tparam InKWeights must be a model of `Range` with value type `InK::FT`.
\tparam HoleIterator must be a model of `InputIterator` with value type being a model of `ConstRange`
                     with value type `SsK::Point_2`.
\tparam HoleWeightsIterator must be a model of `InputIterator` with value type being a model of `ConstRange`
                            with value type `InK::FT`.
\tparam OfKPolygon is a polygon without holes type determined from `OfK` and `InKPolygon`,
                   see Section \ref SLSOffsetPolygonReturnType.

\note If `SsK != OfK` the constructed straight skeleton is converted to `CGAL::Straight_skeleton_2<OfK>`.

\sa `CGAL::create_exterior_skeleton_and_offset_polygons_2()`
\sa `CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2`
*/
template <typename OfKPolygon, typename FT, typename InKPolygon, typename InKWeights,
          typename HoleIterator, typename OfK, typename SsK>
std::vector< std::shared_ptr<OfKPolygon> >
create_interior_weighted_skeleton_and_offset_polygons_2(FT offset,
                                                        const InKPolygon& outer_boundary,
                                                        const InKWeights& outer_boundary_weights,
                                                        HoleIterator holes_begin,
                                                        HoleIterator holes_end,
                                                        HoleWeightsIterator holes_weights_begin,
                                                        HoleWeightsIterator holes_weights_end,
                                                        OfK ofk = CGAL::Exact_predicates_inexact_constructions_kernel,
                                                        SsK ssk = CGAL::Exact_predicates_inexact_constructions_kernel());

/*!
\ingroup PkgStraightSkeleton2WeightedOffsetFunctions

\brief returns a container with all the inner offset polygons at distance `offset` of the 2D polygon `poly`.

Weights of the polygon are given by `weights`, in the same order as the vertices of the polygon:
the `i`-th weight in the range is associated to the contour edge between the `i-1`-th and `i`-th vertices.

\note A temporary straight skeleton is constructed in the interior of the input polygon to obtain the offsets.
The construction of this skeleton is the most expensive operation, therefore, to construct offsets
 at more than one single distance, use the separate functions `create_interior_straight_skeleton_2()`
and `create_offset_polygons_2()` instead.

\tparam OfK must be a model of `Kernel`. It is used to instantiate
            `Polygon_offset_builder_traits_2<OfK>` for constructing the offset polygons.
\tparam SsK must be a model of `Kernel`. It is used to instantiate
            `Straight_skeleton_builder_traits_2<SsK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldNumberType` convertible to `OfK::FT` and `SsK::FT`.
\tparam InKPolygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                   or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).
\tparam InKWeights must be a model of `Range` with value type `InK::FT`.
\tparam OfKPolygon is a polygon without holes type determined by `OfK` and `InKPolygon`,
                   see Section \ref SLSOffsetPolygonReturnType.

\note If `SsK != OfK` the constructed straight skeleton is converted to `CGAL::Straight_skeleton_2<OfK>`.

\sa `CGAL::create_exterior_skeleton_and_offset_polygons_2()`
\sa `CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2`
*/
template <typename OfKPolygon, typename FT, typename InKPolygon, typename InKWeights, typename OfK, typename SsK>
std::vector< std::shared_ptr<OfKPolygon> >
create_interior_weighted_skeleton_and_offset_polygons_2(FT offset,
                                                        const InKPolygon& poly,
                                                        const InKWeights& weights,
                                                        OfK ofk = CGAL::Exact_predicates_inexact_constructions_kernel,
                                                        SsK ssk = CGAL::Exact_predicates_inexact_constructions_kernel());

// ---------------------------------------------- EXTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2WeightedOffsetFunctions

\brief returns a container with all the outer offset polygons at distance `offset` of the 2D polygon `poly`.

Weights of the polygon are given by `weights`, in the same order as the vertices of the polygon:
the `i`-th weight in the range is associated to the contour edge between the `i-1`-th and `i`-th vertices.

The exterior skeleton is limited by an outer rectangular frame placed at a margin sufficiently
large to allow the offset polygons to be constructed.

\note A temporary straight skeleton is constructed in the <I>limited exterior</I> of the input polygon
to obtain the offsets. The construction of this skeleton is the most expensive operation,
therefore, to construct offsets at more than one single distance, use the separate functions
`create_exterior_straight_skeleton_2()` and `create_offset_polygons_2()` instead.

\tparam OfK must be a model of `Kernel`. It is used to instantiate
            `Polygon_offset_builder_traits_2<OfK>` for constructing the offset polygons.
\tparam SsK must be a model of `Kernel`. It is used to instantiate
            `Straight_skeleton_builder_traits_2<SsK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldNumberType` convertible to `OfK::FT` and `SsK::FT`.
\tparam InKPolygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                   or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).
\tparam InKWeights must be a model of `Range` with value type `InK::FT`.
\tparam OfKPolygon is a polygon without holes type determined from `OfK` and `InKPolygon`,
                   see Section \ref SLSOffsetPolygonReturnType.

\note If `SsK != OfK` the constructed straight skeleton is converted to `CGAL::Straight_skeleton_2<OfK>`.

\pre `offset` is positive
\pre poly` is weakly simple, counterclockwise polygon.

\sa `CGAL::create_interior_skeleton_and_offset_polygons_2()`
\sa `CGAL::create_exterior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2`
*/
template <typename OfKPolygon, typename FT, typename InKPolygon, typename InKWeights, typename OfK, typename SsK>
std::vector< std::shared_ptr<OfKPolygon> >
create_exterior_weighted_skeleton_and_offset_polygons_2(FT offset,
                                                        const InKPolygon& poly,
                                                        const InKWeights& weights,
                                                        OfK ofk = Exact_predicates_inexact_constructions_kernel(),
                                                        SsK ssk = Exact_predicates_inexact_constructions_kernel());

} /* namespace CGAL */
