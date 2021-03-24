namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with the offset polygons at distance `offset` obtained from
the straight skeleton `s`.

If `s` is the inner skeleton of a polygon with holes,
the offset polygons will be generated in its interior. If `s` is the outer skeleton
of a polygon with holes, the offset polygons will be generated in its exterior.

\tparam OfK must be a model of `Kernel`. It is used to instantiate
            `Polygon_offset_builder_traits_2<OfK>` for constructing the offset polygons.
\tparam FT must be a model of `FieldNumberType` convertible to `OfK::FT`.
\tparam StraightSkeleton is `Straight_skeleton_2<SsK>`.
\tparam OfKPolygon is a polygon without holes type determined from `OfK`, see Section \ref SLSOffsetPolygonReturnType.

\note If `SsK != OfK` the constructed straight skeleton is converted to `Straight_skeleton_2<OfK>`.

\sa `create_interior_straight_skeleton_2()`
\sa `create_exterior_straight_skeleton_2()`
\sa `create_interior_skeleton_and_offset_polygons_2()`
\sa `create_exterior_skeleton_and_offset_polygons_2()`
\sa `Polygon_offset_builder_2<Ss,Traits,Container>`
*/
template<class OfKPolygon, class FT, class StraightSkeleton, class OfK>
std::vector< boost::shared_ptr<OfKPolygon> >
create_offset_polygons_2 ( FT offset,
                           const StraightSkeleton& s,
                           OfK k = Exact_predicates_inexact_constructions_kernel ) ;

// ---------------------------------------------- EXTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with all the outer offset polygons at distance `offset` of the 2D polygon `poly`.

A temporary straight skeleton is constructed in the <I>limited exterior</I> of the input polygon
to obtain the offsets. The construction of this skeleton is the most
expensive operation, therefore, to construct offsets at more than one
single distance, use the separate functions
`create_exterior_straight_skeleton_2()` and `create_offset_polygons_2()`
instead. The exterior skeleton is limited by an outer rectangular
frame placed at a margin sufficiently large to allow the offset
polygons to be constructed.

\tparam OfK must be a model of `Kernel`. It is used to instantiate
            `Polygon_offset_builder_traits_2<OfK>` for constructing the offset polygons.
\tparam SsK must be a model of `Kernel`. It is used to instantiate
            `Straight_skeleton_builder_traits_2<SsK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldNumberType` convertible to `OfK::FT` and `SsK::FT`.
\tparam InKPolygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`)
                   or a model of `GeneralPolygonWithHoles_2` (e.g. `Polygon_with_holes_2<InK>`).
\tparam OfKPolygon is a polygon without holes type determined from `OfK` and `InKPolygon`,
                   see Section \ref SLSOffsetPolygonReturnType.

\note If `SsK != OfK` the constructed straight skeleton is converted to `Straight_skeleton_2<OfK>`.

\sa `create_interior_skeleton_and_offset_polygons_2()`
\sa `create_exterior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2<Ss,Traits,Container>`
*/
template<class OfKPolygon, class FT, class InKPolygon, class OfK, class SsK>
std::vector< boost::shared_ptr<OfKPolygon> >
create_exterior_skeleton_and_offset_polygons_2( FT offset,
                                                const InKPolygon& poly,
                                                OfK ofk = Exact_predicates_inexact_constructions_kernel,
                                                SsK ssk = Exact_predicates_inexact_constructions_kernel ) ;

// ---------------------------------------------- INTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with all the inner offset polygons at distance `offset` of the 2D polygon
with holes whose outer boundary is `outer_boundary` and its holes are given by `[holes_begin,holes_end)`.

A temporary straight skeleton is constructed in the interior of the input polygon to obtain the offsets.
The construction of this skeleton is the most expensive operation, therefore, to construct offsets
at more than one single distance, use the separate functions `create_interior_straight_skeleton_2()`
and `create_offset_polygons_2()` instead.

\tparam OfK must be a model of `Kernel`. It is used to instantiate
            `Polygon_offset_builder_traits_2<OfK>` for constructing the offset polygons.
\tparam SsK must be a model of `Kernel`. It is used to instantiate
            `Straight_skeleton_builder_traits_2<SsK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldNumberType` convertible to `OfK::FT` and `SsK::FT`.
\tparam HoleIterator must be a model of `InputIterator` with value type being a model of `ConstRange`
                     with value type `SsK::Point_2`.
\tparam InKPolygon must be a model of `SequenceContainer` with value type `InK::Point_2` (e.g. `Polygon_2<InK>`).
\tparam OfKPolygon is a polygon without holes type determined from `OfK` and `InKPolygon`,
                   see Section \ref SLSOffsetPolygonReturnType.

\note If `SsK != OfK` the constructed straight skeleton is converted to `Straight_skeleton_2<OfK>`.

\sa `create_exterior_skeleton_and_offset_polygons_2()`
\sa `create_interior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2<Ss,Traits,Container>`
*/
template<class OfKPolygon, class FT, class InKPolygon, class HoleIterator, class OfK, class SsK>
std::vector< boost::shared_ptr<OfKPolygon> >
create_interior_skeleton_and_offset_polygons_2 ( FT offset,
                                                 const InKPolygon& outer_boundary,
                                                 HoleIterator holes_begin,
                                                 HoleIterator holes_end,
                                                 OfK ofk = CGAL::Exact_predicates_inexact_constructions_kernel,
                                                 SsK ssk = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with all the inner offset polygons at distance `offset` of the 2D polygon `poly`.

A temporary straight skeleton is constructed in the interior of the input polygon to obtain the offsets.
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
\tparam OfKPolygon is a polygon without holes type determined by `OfK` and `InKPolygon`,
                   see Section \ref SLSOffsetPolygonReturnType.

\note If `SsK != OfK` the constructed straight skeleton is converted to `Straight_skeleton_2<OfK>`.

\sa `create_exterior_skeleton_and_offset_polygons_2()`
\sa `create_interior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2<Ss,Traits,Container>`
*/
template<class OfKPolygon, class FT, class InKPolygon, class OfK, class SsK>
std::vector< boost::shared_ptr<OfKPolygon> >
create_interior_skeleton_and_offset_polygons_2 ( FT offset,
                                                 const InKPolygon& poly,
                                                 OfK ofk = CGAL::Exact_predicates_inexact_constructions_kernel,
                                                 SsK ssk = CGAL::Exact_predicates_inexact_constructions_kernel ) ;
} /* namespace CGAL */
