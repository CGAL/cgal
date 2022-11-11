namespace CGAL {

/*!
\ingroup PkgStraightSkeleton2Classes

A straight skeleton event is the simultaneous collision of three offseted oriented straight line segments
`e0*`,`e1*`,`e2*` (`e*` denotes an _offseted_ edge).

This record stores the segments corresponding to the INPUT edges `(e0,e1,e2)` whose offsets intersect
at the event along with their collinearity.

If the event is an edge-event, then `e0*->e1*->e2*` must be consecutive right before the event so that
after the event `e0*` and `e2*` become consecutive. Thus, there are offset vertices `(e0*,e1*)` and `(e1*,e2*)`
in the offset polygon which do not necessarily exist in the original polygon.

If the event is a split-event, `e0*->e1*` must be consecutive right before the event so that after the event
`e0*->right(e2*)` and `left(e2*)->e1*` become consecutive. Thus, there is an offset vertex `(e0*,e1*)` in the
offset polygon which does not necessarily exist in the original polygon.

The offset vertices `(e0*,e1*)` and `(e1*,e2*)` are called the left and right seeds for the event.
A seed is a contour node if the vertex is already present in the input polygon, otherwise is a skeleton node.
If a seed is a skeleton node it is produced by a previous event so it is itself defined as a trisegment, thus,
a trisegment is actually a node in a binary tree.

Since trisegments are tree nodes they must always be handled via the nested smart pointer type `Self_ptr`.

\tparam K_ must be a model of `Kernel`
\tparam Segment_ must be a model of `Kernel::Segment_2`

\note Objects of this type should be constructed using the traits' functor `Construct_ss_trisegment_2`
(see the concept `StraightSkeletonBuilderTraits_2`).

\sa `PolygonOffsetBuilderTraits_2`
*/
template< typename K_, typename Segment_ >
class Trisegment_2
{
public:
  typedef boost::intrusive_ptr<Trisegment_2<K_, Segment_> > Self_ptr ;

  /// returns `e0`.
  const Segment_& e0() const;

  /// returns `e1`.
  const Segment_& e1() const;

  /// returns `e2`.
  const Segment_& e2() const;

}; /* end Trisegment_2 */

} /* end namespace CGAL */
