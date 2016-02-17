namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Functions

returns a container with the offset polygons at distance `offset` obtained from
the straight skeleton `s`. If `s` is the inner skeleton of a polygon with holes,
the offset polygons will be generated in its interior. If `s` is the outer skeleton
of a polygon with holes, the offset polygons will be generated in its exterior.


\tparam K is any \cgal kernel.
\tparam FT is any number type implicitly convertible to `K::FT`.
\tparam StraightSkeleton is `Straight_skeleton_2<K2>`.
If `K != K2` the straight skeleton is converted to `Straight_skeleton_2<K>`.
\tparam Polygon is a model of `VertexContainer_2`.
If this first template parameter is omitted, `Polygon_2` is used.

\sa `create_interior_straight_skeleton_2()`
\sa `create_exterior_straight_skeleton_2()`
\sa `create_interior_skeleton_and_offset_polygons_2()`
\sa `create_exterior_skeleton_and_offset_polygons_2()`
*/

template<class Polygon, class FT, class StraightSkeleton, class K>
std::vector< boost::shared_ptr<Polygon> >
create_offset_polygons_2 ( FT offset
, StraightSkeleton const& s,
, K const& k = Exact_predicates_inexact_constructions_kernel
) ;

/*!
\ingroup PkgStraightSkeleton2Functions

returns a container with all the outer offset polygons at distance
`offset` of the 2D polygon `poly`.  A temporary straight
skeleton is constructed in the `limited exterior` of the input polygon
to obtain the offsets. The construction of this skeleton is the most
expensive operation, therefore, to construct offsets at more than one
single distance, use the separate functions
`create_exterior_straight_skeleton_2()` and `create_offset_polygons_2()`
instead. The exterior skeleton is limited by an outer rectangular
frame placed at a margin sufficiently large to allow the offset
polygons to be constructed.


\tparam OffsettingK must be a \cgal kernel. It is used to instantiate 
    `Polygon_offset_builder_traits_2<K>` for constructing 
    the offset polygons.
\tparam SkeletonK must be a \cgal kernel. It is used to instantiate
   `Straight_skeleton_builder_traits_2<K>` for constructing 
   the straight skeleton.
   If `SkeletonK != OffsettingK` the constructed straight skeleton
   is converted to `Straight_skeleton_2<OffsettingK>`.
\tparam FT must be a number type implicitly convertible to `OffsettingK::FT`.
4. `Straight_skeleton` is `Straight_skeleton_2<K2>`. 
    If `K != K2` the straight skeleton is converted to `Straight_skeleton_2<K>`.

\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`

*/
template<class FT, class Polygon, class OffsettingK, class SkeletonK>
std::vector< boost::shared_ptr<Polygon> >
create_exterior_skeleton_and_offset_polygons_2
  ( FT             offset
  , Polygon const& poly
  , OffsettingK    ofk  = Exact_predicates_inexact_constructions_kernel
  , SkeletonK      ssk  = Exact_predicates_inexact_constructions_kernel
  ) ;

/*!
\ingroup PkgStraightSkeleton2Functions
returns a container with all the inner offset polygons at distance `offset` of the 2D polygon with holes whose outer boundary is `outer_boundary` and its holes are given by `[holes_begin,holes_end]`.
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
\tparam Straight_skeleton is `Straight_skeleton_2<K2>`. 
   If `K != K2` the straight skeleton is converted to `Straight_skeleton_2<K>`.
\tparam HoleIterator is an `InputIterator` with `value_type` being
        `CGAL::Polygon_2<OffsettingK>` or a standard container of
        `OffsettingK::Point_2` elements 
\tparam Polygon is an `InputIterator` with `value_type` being
        `CGAL::Polygon_2<OffsettingK>` or a standard container of
        `OffsettingK::Point_2` elements 

\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2`

*/
template<class FT, class Polygon, class HoleIterator, class OffsettingK, class SkeletonK>
std::vector< boost::shared_ptr<Polygon> >
create_interior_skeleton_and_offset_polygons_2 ( FT             offset
                                               , Polygon const& outer_boundary
                                               , HoleIterator   holes_begin
                                               , HoleIterator   holes_end
                                               , OffsettingK   ofk 
                                                  = CGAL::Exact_predicates_inexact_constructions_kernel
                                               , SkeletonK      ssk 
                                                  = CGAL::Exact_predicates_inexact_constructions_kernel
                                               ) ;

/*!
\ingroup PkgStraightSkeleton2Functions

returns a container with all the inner offset polygons at distance `offset` of the 2D polygon `poly`.
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
\tparam Polygon is an `InputIterator` with `value_type` being
        `CGAL::Polygon_2<OffsettingK>` or a standard container of
        `OffsettingK::Point_2` elements 


\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2` 
*/
template<class FT, class Polygon, class OffsettingK, class SkeletonK>
std::vector< boost::shared_ptr<Polygon> >
create_interior_skeleton_and_offset_polygons_2 ( FT             offset
                                               , Polygon const& poly
                                               , OffsettingK     ofk 
                                                  = CGAL::Exact_predicates_inexact_constructions_kernel
                                               , SkeletonK      ssk 
                                                  = CGAL::Exact_predicates_inexact_constructions_kernel
                                               ) ;
} /* namespace CGAL */
