namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Functions

returns a new `Straight_skeleton_2<K>` in the interior of the 2D polygon with holes `poly_with_holes`.

\sa `create_exterior_straight_skeleton_2()`
\sa `Straight_skeleton_builder_2` 

*/
template<class K, class C>
boost::shared_ptr< Straight_skeleton_2<K> >
  create_interior_straight_skeleton_2 ( Polygon_with_holes<K,C> poly_with_holes
                                    , K const&  k = Exact_predicates_inexact_constructions_kernel
                                    ) ;
} /* namespace CGAL */
