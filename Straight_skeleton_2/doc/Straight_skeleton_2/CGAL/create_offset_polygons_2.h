namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with the offset polygons at distance `offset` obtained from
the straight skeleton `s`.

If `s` is the inner skeleton of a polygon with holes,
the offset polygons will be generated in its interior. If `s` is the outer skeleton
of a polygon with holes, the offset polygons will be generated in its exterior.

\tparam OffsettingK must be a model of `Kernel`. It is used to instantiate
                    `Polygon_offset_builder_traits_2<OffsettingK>` for constructing the offset polygons.
\tparam FT must be a model of `FieldType` convertible to `OffsettingK::FT`.
\tparam StraightSkeleton is `Straight_skeleton_2<SkeletonK>`.
\tparam Polygon must be a model of `SequenceContainer` with value type `OffsettingK::Point_2`
                (e.g. `Polygon_2<OffsettingK>` or `Polygon_with_holes_2<OffsettingK>`).

\note If `SkeletonK != OffsettingK` the constructed straight skeleton is converted to `Straight_skeleton_2<OffsettingK>`.

\sa `create_interior_straight_skeleton_2()`
\sa `create_exterior_straight_skeleton_2()`
\sa `create_interior_skeleton_and_offset_polygons_2()`
\sa `create_exterior_skeleton_and_offset_polygons_2()`
\sa `Polygon_offset_builder_2`
*/
template<class Polygon, class FT, class StraightSkeleton, class OffsettingK>
std::vector< boost::shared_ptr<Polygon> >
create_offset_polygons_2 ( FT offset,
                           const StraightSkeleton& s,
                           OffsettingK k = Exact_predicates_inexact_constructions_kernel ) ;

// ---------------------------------------------- EXTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with all the outer offset polygons at distance `offset` of the 2D polygon `poly`.

A temporary straight skeleton is constructed in the `limited exterior` of the input polygon
to obtain the offsets. The construction of this skeleton is the most
expensive operation, therefore, to construct offsets at more than one
single distance, use the separate functions
`create_exterior_straight_skeleton_2()` and `create_offset_polygons_2()`
instead. The exterior skeleton is limited by an outer rectangular
frame placed at a margin sufficiently large to allow the offset
polygons to be constructed.

\tparam OffsettingK must be a model of `Kernel`. It is used to instantiate
                    `Polygon_offset_builder_traits_2<OffsettingK>` for constructing the offset polygons.
\tparam SkeletonK must be a model of `Kernel`. It is used to instantiate
                  `Straight_skeleton_builder_traits_2<SkeletonK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldType` convertible to `OffsettingK::FT` and `SkeletonK::FT`.
\tparam Polygon must be a model of `SequenceContainer` with value type `OffsettingK::Point_2`
                (e.g. `Polygon_2<OffsettingK>` or `Polygon_with_holes_2<OffsettingK>`).

\note If `SkeletonK != OffsettingK` the constructed straight skeleton is converted to `Straight_skeleton_2<OffsettingK>`.

\sa `create_interior_skeleton_and_offset_polygons_2()`
\sa `create_exterior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2`
*/
template<class FT, class Polygon, class OffsettingK, class SkeletonK>
std::vector< boost::shared_ptr<Polygon> >
create_exterior_skeleton_and_offset_polygons_2( FT offset,
                                                const Polygon& poly,
                                                OffsettingK ofk = Exact_predicates_inexact_constructions_kernel, // @fixme EPICK?
                                                SkeletonK ssk = Exact_predicates_inexact_constructions_kernel ) ;

// ---------------------------------------------- INTERIOR -----------------------------------------

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with all the inner offset polygons at distance `offset` of the 2D polygon
with holes whose outer boundary is `outer_boundary` and its holes are given by `[holes_begin,holes_end]`.

A temporary straight skeleton is constructed in the interior of the input polygon to obtain the offsets.
The construction of this skeleton is the most expensive operation, therefore, to construct offsets
at more than one single distance, use the separate functions `create_interior_straight_skeleton_2()`
and `create_offset_polygons_2()` instead.

\tparam OffsettingK must be a model of `Kernel`. It is used to instantiate
                    `Polygon_offset_builder_traits_2<OffsettingK>` for constructing the offset polygons.
\tparam SkeletonK must be a model of `Kernel`. It is used to instantiate
                  `Straight_skeleton_builder_traits_2<SkeletonK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldType` convertible to `OffsettingK::FT` and `SkeletonK::FT`.
\tparam HoleIterator must be a model of `InputIterator` with value type being a model of `ConstRange`
                     with value type `SkeletonK::Point_2`.
\tparam Polygon must be a model of `SequenceContainer` with value type `SkeletonK::Point_2`
                (e.g. `Polygon_2<SkeletonK>` or `Polygon_with_holes_2<SkeletonK>`).

\note If `SkeletonK != OffsettingK` the constructed straight skeleton is converted to `Straight_skeleton_2<OffsettingK>`.

\sa `create_exterior_skeleton_and_offset_polygons_2()`
\sa `create_interior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2`
*/
template<class FT, class Polygon, class HoleIterator, class OffsettingK, class SkeletonK>
std::vector< boost::shared_ptr<Polygon> >
create_interior_skeleton_and_offset_polygons_2 ( FT offset,
                                                 const Polygon& outer_boundary,
                                                 HoleIterator holes_begin,
                                                 HoleIterator holes_end,
                                                 OffsettingK ofk = CGAL::Exact_predicates_inexact_constructions_kernel,
                                                 SkeletonK ssk = CGAL::Exact_predicates_inexact_constructions_kernel ) ;

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

returns a container with all the inner offset polygons at distance `offset` of the 2D polygon `poly`.

A temporary straight skeleton is constructed in the interior of the input polygon to obtain the offsets.
The construction of this skeleton is the most expensive operation, therefore, to construct offsets
 at more than one single distance, use the separate functions `create_interior_straight_skeleton_2()`
and `create_offset_polygons_2()` instead.

\tparam OffsettingK is the \cgal kernel used to instantiate
                    `Polygon_offset_builder_traits_2<OffsettingK>` for constructing the offset polygons.
\tparam SkeletonK is the \cgal kernel used to instantiate
                  `Straight_skeleton_builder_traits_2<SkeletonK>` for constructing the straight skeleton.
\tparam FT must be a model of `FieldType` convertible to `OffsettingK::FT` and `SkeletonK::FT`.
\tparam Polygon must be a model of `SequenceContainer` with value type `SkeletonK::Point_2`
                (e.g. `Polygon_2<SkeletonK>` or `Polygon_with_holes_2<SkeletonK>`).

\note If `SkeletonK != OffsettingK` the constructed straight skeleton is converted to `Straight_skeleton_2<OffsettingK>`.

\sa `create_exterior_skeleton_and_offset_polygons_2()`
\sa `create_interior_skeleton_and_offset_polygons_with_holes_2()`
\sa `Polygon_offset_builder_2`
*/
template<class FT, class Polygon, class OffsettingK, class SkeletonK>
std::vector< boost::shared_ptr<Polygon> >
create_interior_skeleton_and_offset_polygons_2 ( FT offset,
                                                 const Polygon& poly,
                                                 OffsettingK ofk = CGAL::Exact_predicates_inexact_constructions_kernel,
                                                 SkeletonK ssk = CGAL::Exact_predicates_inexact_constructions_kernel ) ;
} /* namespace CGAL */
