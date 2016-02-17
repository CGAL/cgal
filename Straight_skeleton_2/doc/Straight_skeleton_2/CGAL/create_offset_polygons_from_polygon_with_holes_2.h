namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Functions

returns a container with all the outer offset polygons <I>with
holes</I> at distance `offset` of the 2D polygon `poly`.

This is equivalent to `arrange_offset_polygons_2(create_exterior_skeleton_and_offset_polygons_2(offset,poly,okk,ssk))`.
*/
template<class FT, class OffsettingK, class SkeletonK, class C>
std::vector< boost::shared_ptr< Polygon_with_holes_2<OffsettingK,C> > >
create_exterior_skeleton_and_offset_polygons_with_holes_2 
  ( FT offset
  , CGAL::Polygon_2<OffsettingK,C> poly_with_holes
  , OffsettingK ofk = Exact_predicates_inexact_constructions_kernel
  , SkeletonK ssk   = Exact_predicates_inexact_constructions_kernel
  ) ;

/*!
\ingroup PkgStraightSkeleton2Functions

returns a container with all the inner offset polygons at distance `offset` of the 2D polygon with holes `poly_with_holes`.
A temporary straight skeleton is constructed in the interior of the input polygon to obtain the offsets. The construction of this skeleton is the most expensive operation, therefore, to construct offsets at more than one single distance, use the separate functions `create_interior_straight_skeleton_2()` and `create_offset_polygons_2()` instead.

\tparam OffsettingK is the \cgal kernel used to instantiate
   `Polygon_offset_builder_traits_2<K>` for constructing 
   the offset polygons.
\tparam SkeletonK is the \cgal kernel used to instantiate
   `Straight_skeleton_builder_traits_2<K>` for constructing 
   the straight skeleton.
   If `SkeletonK != OffsettingK` the constructed straight skeleton
   is converted to `Straight_skeleton_2<OffsettingK>`.
\tparam FT is any number type implicitly convertible to `OffsettingK::FT`.
\tparam C is the container of points in the polygons.


\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`

*/
template<class FT, class OffsettingK, class SkeletonK, class C>
std::vector< boost::shared_ptr< Polygon_2<OffsettingK,C> > >
create_interior_skeleton_and_offset_polygons_2 ( FT offset
                                               , CGAL::Polygon_with_holes<OffsettingK,C> poly_with_holes
                                               , OffsettingK ofk 
                                                  = CGAL::Exact_predicates_inexact_constructions_kernel
                                               , SkeletonK ssk 
                                                  = CGAL::Exact_predicates_inexact_constructions_kernel
                                    ) ;

/*!
\ingroup PkgStraightSkeleton2Functions

returns a container with all the inner offset polygons <I>with holes</I> at distance `offset` of the 2D polygon `poly`.
This is equivalent to `arrange_offset_polygons_2(create_interior_skeleton_and_offset_polygons_2(offset,poly,okk,ssk))`.

\tparam OffsettingK is the \cgal kernel used to instantiate
   `Polygon_offset_builder_traits_2<K>` for constructing 
   the offset polygons.
\tparam SkeletonK is the \cgal kernel used to instantiate
   `Straight_skeleton_builder_traits_2<K>` for constructing 
   the straight skeleton.
   If `SkeletonK != OffsettingK` the constructed straight skeleton
   is converted to `Straight_skeleton_2<OffsettingK>`.
\tparam FT is any number type implicitly convertible to `OffsettingK::FT`.
\tparam C is the container of points in the polygons.


\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`

*/
template<class FT, class OffsettingK, class SkeletonK, class C>
std::vector< boost::shared_ptr< Polygon_with_holes_2<OffsettingK,C> > >
create_interior_skeleton_and_offset_polygons_with_holes_2 
  ( FT offset
  , CGAL::Polygon_2<OffsettingK,C> poly
  , OffsettingK ofk = CGAL::Exact_predicates_inexact_constructions_kernel
  , SkeletonK ssk   = CGAL::Exact_predicates_inexact_constructions_kernel
  ) ;

/*!
\ingroup PkgStraightSkeleton2Functions

returns a container with all the inner offset polygons <I>with holes</I> at distance `offset` of the 2D polygon with holes `poly_with_holes`.
This is equivalent to `arrange_offset_polygons_2(create_interior_skeleton_and_offset_polygons_2(offset,poly_with_holes,okk,ssk))`.

\tparam OffsettingK is the \cgal kernel used to instantiate
   `Polygon_offset_builder_traits_2<K>` for constructing 
   the offset polygons.
\tparam SkeletonK is the \cgal kernel used to instantiate
   `Straight_skeleton_builder_traits_2<K>` for constructing 
   the straight skeleton.
   If `SkeletonK != OffsettingK` the constructed straight skeleton
   is converted to `Straight_skeleton_2<OffsettingK>`.
\tparam FT is any number type implicitly convertible to `OffsettingK::FT`.
\tparam C is the container of points in the polygons.
\tparam Straight_skeleton is `Straight_skeleton_2<K2>`. 
   If `K != K2` the straight skeleton is converted to `Straight_skeleton_2<K>`.

\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`

*/
template<class FT, class OffsettingK, class SkeletonK, class C>
std::vector< boost::shared_ptr< Polygon_with_holes_2<OffsettingK,C> > >
create_interior_skeleton_and_offset_polygons_with_holes_2 
  ( FT offset
  , CGAL::Polygon_with_holes_2<OffsettingK,C> poly_with_holes
  , OffsettingK ofk = CGAL::Exact_predicates_inexact_constructions_kernel
  , SkeletonK ssk   = CGAL::Exact_predicates_inexact_constructions_kernel
  ) ;

} /* namespace CGAL */
