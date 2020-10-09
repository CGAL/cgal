namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2OffsetFunctions

The function `arrange_offset_polygons_2()` arranges the sequence of 2D polygons obtained
by `create_offset_polygons_2()` into 2D polygons with holes by determining geometric parent-hole
relationships using a simple algorithm based on the particular characteristics of offset polygons.

The function determines parent-hole relationships among the polygons given by `[begin,end]` creating
`boost::shared_ptr< GeneralPolygonWithHoles_2 >` objects added to the output sequence given `out`.
A `CLOCKWISE` oriented polygon `H` is a hole of a `COUNTERCLOCKWISE` polygon `P`, iff at least one vertex of `H` is `ON_BOUNDED_SIDE` of `P`.

This function should not be used to arrange arbitrary polygons into polygons with holes unless they meet the requirements specified below.

\tparam K must be a model of `Kernel`.
\tparam InputPolygonPtrIterator must be a model of `InputIterator` whose `value_type` is a smart pointer
(such as `boost::shared_ptr`) whose `element_type` is a model of `SequenceContainer` with value type `K::Point_2`.
\tparam OutputPolygonWithHolesPtrIterator must be a model of `OutputIterator` whose `value_type` is a smart pointer
(such as `boost::shared_ptr`) whose `element_type` is a model of `GeneralPolygonWithHoles_2`.

\pre The input polygons must be simple.
\pre The set of input polygons are unique and interior disjoint. That is, given distinct polygons
`P` and `Q`, there are vertices of `P` which are not on the boundary of `Q` and are all on the
bounded or unbounded side of `Q` (but not both).

\return `true` if no error was encountered, and `false` otherwise.

\sa `create_exterior_straight_skeleton_2()`
\sa `CGAL::Straight_skeleton_builder_2<Traits,Ss,Visitor>`

*/
template<class K, class InputPolygonPtrIterator, class OutputPolygonWithHolesPtrIterator>
bool arrange_offset_polygons_2 ( InputPolygonPtrIterator begin,
                                 InputPolygonPtrIterator end,
                                 OutputPolygonWithHolesPtrIterator out,
                                 const K& k) ;

} /* namespace CGAL */
