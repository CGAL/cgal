namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with all the inner offset polygons <I>with holes</I> at distance `offset`
of the 2D polygon with holes `poly_with_holes`.

This is equivalent to `arrange_offset_polygons_2(create_interior_skeleton_and_offset_polygons_2(offset, poly_with_holes, ofk, ssk))`.

\tparam OffsettingK must be a model of `Kernel`. It is used to instantiate
                    `Polygon_offset_builder_traits_2<OffsettingK>` for constructing the offset polygons.
\tparam SkeletonK must be a model of `Kernel`. It is used to instantiate
                  `Straight_skeleton_builder_traits_2<SkeletonK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldType` convertible to `OffsettingK::FT` and `SkeletonK::FT`.
\tparam Polygon must be `Polygon_2<OffsettingK>` or `Polygon_with_holes_2<OffsettingK>`.

\note If `SkeletonK != OffsettingK` the constructed straight skeleton is converted to `Straight_skeleton_2<OffsettingK>`.

\sa `create_interior_skeleton_and_offset_polygons_2()`
\sa `create_exterior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2`
*/
template<class FT, class Polygon, class OffsettingK, class SkeletonK>
std::vector< boost::shared_ptr< Polygon > >
create_interior_skeleton_and_offset_polygons_with_holes_2(FT offset,
                                                          const Polygon& poly_with_holes,
                                                          OffsettingK ofk = CGAL::Exact_predicates_inexact_constructions_kernel,
                                                          SkeletonK ssk = CGAL::Exact_predicates_inexact_constructions_kernel);


// ---------------------------------------------- EXTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with all the outer offset polygons <I>with holes</I>
at distance `offset` of the 2D polygon `poly_with_holes`.

This is equivalent to `arrange_offset_polygons_2(create_exterior_skeleton_and_offset_polygons_2(offset, poly_with_holes, ofk, ssk))`.

\tparam OffsettingK must be a model of `Kernel`. It is used to instantiate
                    `Polygon_offset_builder_traits_2<OffsettingK>` for constructing the offset polygons.
\tparam SkeletonK must be a model of `Kernel`. It is used to instantiate
                  `Straight_skeleton_builder_traits_2<SkeletonK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldType` convertible to `OffsettingK::FT` and `SkeletonK::FT`.
\tparam Polygon must be `Polygon_2<OffsettingK>` or `Polygon_with_holes_2<OffsettingK>`.

\note If `SkeletonK != OffsettingK` the constructed straight skeleton is converted to `Straight_skeleton_2<OffsettingK>`.

\sa `create_exterior_skeleton_and_offset_polygons_2()`
\sa `create_interior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2`
*/
template<class FT,  class Polygon, class OffsettingK, class SkeletonK>
std::vector<boost::shared_ptr<Polygon> >
create_exterior_skeleton_and_offset_polygons_with_holes_2(FT offset,
                                                          const Polygon& poly_with_holes,
                                                          OffsettingK ofk = Exact_predicates_inexact_constructions_kernel,
                                                          SkeletonK ssk = Exact_predicates_inexact_constructions_kernel);

} /* namespace CGAL */
