namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with all the inner offset polygons <I>with holes</I> at distance `offset`
of the 2D polygon with holes `poly_with_holes`.

This is equivalent to `arrange_offset_polygons_2(create_interior_skeleton_and_offset_polygons_2(offset, poly_with_holes, ofk, ssk))`.

\tparam OfK must be a model of `Kernel`. It is used to instantiate
            `Polygon_offset_builder_traits_2<OfK>` for constructing the offset polygons.
\tparam SsK must be a model of `Kernel`. It is used to instantiate
            `Straight_skeleton_builder_traits_2<SsK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldNumberType` convertible to `OfK::FT` and `SsK::FT`.
\tparam InKPolygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                   or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).
\tparam OfKPolygon is a polygon without holes type determined by `OfK` and `InKPolygon`,
                   see Section \ref SLSOffsetPolygonReturnType.

\note If `SsK != OfK` the constructed straight skeleton is converted to `Straight_skeleton_2<OfK>`.

\sa `create_interior_skeleton_and_offset_polygons_2()`
\sa `create_exterior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2<Ss,Traits,Container>`
*/
template<class OfKPolygon, class FT, class InKPolygon, class OfK, class SsK>
std::vector< boost::shared_ptr< OfKPolygon > >
create_interior_skeleton_and_offset_polygons_with_holes_2(FT offset,
                                                          const InKPolygon& poly_with_holes,
                                                          OfK ofk = CGAL::Exact_predicates_inexact_constructions_kernel,
                                                          SsK ssk = CGAL::Exact_predicates_inexact_constructions_kernel);


// ---------------------------------------------- EXTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with all the outer offset polygons <I>with holes</I>
at distance `offset` of the 2D polygon `poly_with_holes`.

This is equivalent to `arrange_offset_polygons_2(create_exterior_skeleton_and_offset_polygons_2(offset, poly_with_holes, ofk, ssk))`.

\tparam OfK must be a model of `Kernel`. It is used to instantiate
            `Polygon_offset_builder_traits_2<OfK>` for constructing the offset polygons.
\tparam SsK must be a model of `Kernel`. It is used to instantiate
            `Straight_skeleton_builder_traits_2<SsK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldNumberType` convertible to `OfK::FT` and `SsK::FT`.
\tparam InKPolygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                   or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).
\tparam OfKPolygon is a polygon without holes type determined by `OfK` and `InKPolygon`,
                   see Section \ref SLSOffsetPolygonReturnType.

\note If `SsK != OfK` the constructed straight skeleton is converted to `Straight_skeleton_2<OfK>`.

\sa `create_exterior_skeleton_and_offset_polygons_2()`
\sa `create_interior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2<Ss,Traits,Container>`
*/
template<class OfKPolygon, class FT, class InKPolygon, class OfK, class SsK>
std::vector<boost::shared_ptr<OfKPolygon> >
create_exterior_skeleton_and_offset_polygons_with_holes_2(FT offset,
                                                          const InKPolygon& poly_with_holes,
                                                          OfK ofk = Exact_predicates_inexact_constructions_kernel,
                                                          SsK ssk = Exact_predicates_inexact_constructions_kernel);

} /* namespace CGAL */
