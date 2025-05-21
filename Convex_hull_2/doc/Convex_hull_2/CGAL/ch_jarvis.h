namespace CGAL {

/*!
\ingroup PkgConvexHull2Functions

generates the counterclockwise sequence of extreme points
of the points in the range [`first`,`beyond`).
The resulting sequence is placed starting at position
`result`, and the past-the-end iterator for the resulting
sequence is returned. It is not specified at which point the
cyclic sequence of extreme points is cut into a linear sequence.

\pre The source range [`first`,`beyond`) does not contain `result`.

The default traits class `Default_traits` is the kernel in which the
value type of `ForwardIterator` is defined.

\cgalHeading{Requirements}

<OL>
<LI>The value type of `ForwardIterator` and
`OutputIterator` is equivalent to `Traits::Point_2`.
<LI>`Traits` defines the following subset of types from
the concept `ConvexHullTraits_2` and their corresponding member
functions that return instances of these types:
<UL>
<LI>`Traits::Point_2`,
<LI>`Traits::Equal_2`,
<LI>`Traits::Less_rotate_ccw_2`,
<LI>`Traits::Less_xy_2`.
</UL>
</OL>

\sa `CGAL::ch_akl_toussaint()`
\sa `CGAL::ch_bykat()`
\sa `CGAL::ch_eddy()`
\sa `CGAL::ch_graham_andrew()`
\sa `CGAL::ch_jarvis_march()`
\sa `CGAL::ch_melkman()`
\sa `CGAL::lower_hull_points_2()`
\sa `CGAL::upper_hull_points_2()`
\sa `CGAL::convex_hull_2()`

\cgalHeading{Implementation}

This function uses the Jarvis march (gift-wrapping)
algorithm \cgalCite{j-ichfs-73}. This algorithm requires \cgalBigO{n h} time
in the worst case for \f$ n\f$ input points with \f$ h\f$ extreme points.

*/
template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
ch_jarvis(ForwardIterator first, ForwardIterator beyond,
          OutputIterator result,
          const Traits& ch_traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgConvexHull2Subsequence

generates the counterclockwise sequence of extreme
points from a given set of input points that line between two input points.

More precisely, it generates the counterclockwise subsequence of
extreme points between `start_p` and `stop_p` of the
points in the range [`first`,`beyond`), starting at
position `result` with point `start_p`. The last point
generated is the point preceding `stop_p` in the
counterclockwise order of extreme points.

The default traits class `Default_traits` is the kernel in which the
type `ForwardIterator::value_type` is defined.

\cgalHeading{Requirements}

<OL>
<LI>`ForwardIterator::value_type` and
`OutputIterator::value_type`
are equivalent to `Traits::Point_2`.
<LI>`Traits` defines the following subset of types from
the concept `ConvexHullTraits_2` and their corresponding member
functions that return instances of these types:
<UL>
<LI>`Traits::Point_2`,
<LI>`Traits::Equal_2`,
<LI>`Traits::Less_rotate_ccw_2`.
</UL>
</OL>

\sa `CGAL::ch_jarvis()`
\sa `CGAL::lower_hull_points_2()`
\sa `CGAL::upper_hull_points_2()`

\cgalHeading{Implementation}

The function uses the Jarvis march (gift-wrapping) algorithm \cgalCite{j-ichfs-73}.
This algorithm requires \cgalBigO{n h} time in the worst
case for \f$ n\f$ input points with \f$ h\f$ extreme points

\pre `start_p` and `stop_p` are extreme points with respect to
the points in the range [`first`,`beyond`) and `stop_p`
is an element of range [`first`,`beyond`).
*/
template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
ch_jarvis_march(ForwardIterator first, ForwardIterator beyond,
                const Traits::Point_2& start_p,
                const Traits::Point_2& stop_p,
                OutputIterator result,
                const Traits& ch_traits = Default_traits);

} /* namespace CGAL */
