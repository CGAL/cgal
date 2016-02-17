namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Functions

Computes the separation required between a polygon and the outer frame used to obtain an exterior skeleton suitable for the computation of outer offset polygons at a given distance. 

Given a non-degenerate strictly-simple 2D polygon whose vertices are passed 
in the range [`first`,`beyond`), calculates the largest euclidean distance
`d` between each input vertex and its corresponding offset vertex at
a distance `offset`.

If such a distance can be approximately computed, returns an `optional<FT>` with the value `d + (offset * 1.05)`. If the distance cannot be computed, not even approximately, due to overflow for instance, returns an empty `optional<FT>` (an <I>absent result</I>).

This result is the required separation between the input polygon
and the rectangular frame used to construct an exterior offset contour
at distance `offset` (which is done by placing the polygon as a hole of that frame).

Such a separation must be computed in this way because if the frame is 
too close to the polygon, the inward offset contour from the frame could
collide with the outward offset contour of the polygon, resulting in a merged
contour offset instead of two contour offsets, one of them corresponding to the frame.

Simply using `2*offset` as the separation is incorrect since `offset` is the distance 
between an offset line and its original, not between an offset vertex and its original.
The later, which is calculated by this function and needed to place the frame sufficiently
away from the polygon, can be thousands of times larger than `offset`.

If the result is <I>absent</I>, any attempt to construct an exterior offset polygon at distance `offset` will fail. This will occur whenever the polygon has a vertex with an internal angle approaching `0` (because the offset vertex of a vertex whose internal angle equals 0 is at <I>infinity</I> ).

\pre `offset > 0`.
\pre The range [`first`,`beyond`) contains the vertices of a non-degenerate strictly-simple 2D polygon.

The default traits class `Default_traits` is an instance of the 
class `Polygon_offset_builder_traits_2<Kernel>` parameterized on 
the kernel in which the type `InputIterator::value_type` is defined. 


\tparam InputIterator must have a `value_type` equivalent to `Traits::Point_2`. 
\tparam Traits must be a model for `PolygonOffsetBuilderTraits_2`.

\sa `PolygonOffsetBuilderTraits_2` 
\sa `Polygon_offset_builder_traits_2` 


*/
template <class InputIterator, class Traits>
boost::optional< typename Traits::FT >
compute_outer_frame_margin( InputIterator first
, InputIterator beyond
, typename Traits::FT offset
, Traits const& traits = Default_traits
);

} /* namespace CGAL */
