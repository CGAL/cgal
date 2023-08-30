namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

\brief returns a container with all the inner offset polygons <I>with holes</I> at distance `offset`
of the 2D polygon with holes `poly_with_holes`.

This is equivalent to `arrange_offset_polygons_2(create_interior_weighted_skeleton_and_offset_polygons_2(offset, poly_with_holes, ofk, ssk))`.

\tparam OfK must be a model of `Kernel`. It is used to instantiate
            `Polygon_offset_builder_traits_2<OfK>` for constructing the offset polygons.
\tparam SsK must be a model of `Kernel`. It is used to instantiate
            `Straight_skeleton_builder_traits_2<SsK>` for constructing the weighted straight skeleton.
\tparam FT must be a model of `FieldNumberType` convertible to `OfK::FT` and `SsK::FT`.
\tparam InKPolygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                   or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).
\tparam OfKPolygon is a polygon without holes type determined by `OfK` and `InKPolygon`,
                   see Section \ref SLSOffsetPolygonReturnType.

\note If `SsK != OfK` the constructed weighted straight skeleton is converted to `CGAL::Straight_skeleton_2<OfK>`.

\sa `CGAL::create_interior_weighted_skeleton_and_offset_polygons_2()`
\sa `CGAL::create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2`
*/
template<class OfKPolygon, class FT, class InKPolygon, class OfK, class SsK>
std::vector< std::shared_ptr< OfKPolygon > >
create_interior_weighted_skeleton_and_offset_polygons_with_holes_2(FT offset,
                                                                   const InKPolygon& poly_with_holes,
                                                                   OfK ofk = CGAL::Exact_predicates_inexact_constructions_kernel,
                                                                   SsK ssk = CGAL::Exact_predicates_inexact_constructions_kernel);


// ---------------------------------------------- EXTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

\brief returns a container with all the outer offset polygons <I>with holes</I>
at distance `offset` of the 2D polygon `poly_with_holes`. Note that the
offset of the outer frame is ignored.

This is equivalent to a call to `CGAL::arrange_offset_polygons_2()` on the
output of \link CGAL::create_exterior_weighted_skeleton_and_offset_polygons_2() `create_exterior_weighted_skeleton_and_offset_polygons_2(offset, poly_with_holes, ofk, ssk))` \endlink
after having filtered out the polygon corresponding to the offset of the outer frame and
having reversed the orientation of all other polygons.

\tparam OfK must be a model of `Kernel`. It is used to instantiate
            `Polygon_offset_builder_traits_2<OfK>` for constructing the offset polygons.
\tparam SsK must be a model of `Kernel`. It is used to instantiate
            `Straight_skeleton_builder_traits_2<SsK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldNumberType` convertible to `OfK::FT` and `SsK::FT`.
\tparam InKPolygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                   or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).
\tparam OfKPolygon is a polygon without holes type determined by `OfK` and `InKPolygon`,
                   see Section \ref SLSOffsetPolygonReturnType.

\note If `SsK != OfK` the constructed weighted straight skeleton is converted to `CGAL::Straight_skeleton_2<OfK>`.

\sa `CGAL::create_exterior_weighted_skeleton_and_offset_polygons_2()`
\sa `CGAL::create_interior_weighted_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2`
*/
template<class OfKPolygon, class FT, class InKPolygon, class OfK, class SsK>
std::vector<std::shared_ptr<OfKPolygon> >
create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2(FT offset,
                                                                   const InKPolygon& poly_with_holes,
                                                                   OfK ofk = Exact_predicates_inexact_constructions_kernel(),
                                                                   SsK ssk = Exact_predicates_inexact_constructions_kernel());

} /* namespace CGAL */
