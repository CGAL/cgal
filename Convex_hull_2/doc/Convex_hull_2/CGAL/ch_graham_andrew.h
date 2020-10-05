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
value type of  `InputIteratore` is defined.

\cgalHeading{Requirements}

<OL>
<LI>The value type of `InputIterator` and
`OutputIterator` is equivalent to `Traits::Point_2`.
<LI>`Traits` defines the following subset of types from
the concept `ConvexHullTraits_2` and their corresponding member
functions that return instances of these types:
<UL>
<LI>`Traits::Point_2`,
<LI>`Traits::Less_xy_2`,
<LI>`Traits::Left_turn_2`,
<LI>`Traits::Equal_2`.
</UL>
</OL>

\sa `CGAL::ch_akl_toussaint()`
\sa `CGAL::ch_bykat()`
\sa `CGAL::ch_eddy()`
\sa `CGAL::ch_graham_andrew_scan()`
\sa `CGAL::ch_jarvis()`
\sa `CGAL::ch_melkman()`
\sa `CGAL::convex_hull_2()`
\sa `CGAL::lower_hull_points_2()`
\sa `CGAL::upper_hull_points_2()`

\cgalHeading{Implementation}

This function implements Andrew's variant of the Graham
scan algorithm \cgalCite{a-aeach-79} and follows the presentation of Mehlhorn
\cgalCite{m-mdscg-84}. This algorithm requires \f$ O(n \log n)\f$ time
in the worst case for \f$ n\f$ input points.


*/
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_graham_andrew( InputIterator first,
InputIterator beyond,
OutputIterator result,
const Traits & ch_traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgConvexHull2Subsequence

generates the counterclockwise sequence of extreme
points from a given sequence of input points that are not left of the line defined
by the first and last points in this sequence.

More precisely, it generates the counterclockwise sequence of extreme
points from a given sequence of input points that are not left of the line
\f$ pq\f$ defined by the first (\f$p\f$) and last (\f$q\f$) points
in this sequence (\f$ p\f$ is the value of `first` and \f$ q\f$ is
the value of `beyond` \f$ -1\f$).
The resulting sequence is placed
starting at `result` with \f$ p\f$; point \f$ q\f$ is omitted.
The past-the-end iterator for the sequence is returned.

\pre The range [`first`,`beyond`) contains at least two different points. The points in [`first`,`beyond`) are sorted with respect to \f$ pq\f$, <I>i.e.</I>, the sequence of points in [`first`,`beyond`) define a counterclockwise polygon, for which the Graham-Sklansky-procedure \cgalCite{s-mcrm-72} works.

The default traits class `Default_traits` is the kernel in which the
type `BidirectionalIterator::value_type` is defined.

\cgalHeading{Requirements}

<OL>
<LI>`BidirectionalIterator::value_type` and
`OutputIterator::value_type`
are equivalent to `Traits::Point_2`.
<LI>`Traits` defines the following two types from
the concept `ConvexHullTraits_2` and their corresponding member
functions that return instances of these types:
<UL>
<LI>`Traits::Point_2`,
<LI>`Traits::Left_turn_2`.
</UL>
</OL>

\sa `CGAL::ch_graham_andrew()`
\sa `CGAL::lower_hull_points_2()`
\sa `CGAL::upper_hull_points_2()`

\cgalHeading{Implementation}

This algorithm requires \f$ O(n)\f$ time in the worst case for
\f$ n\f$ input points.

\cgalHeading{Example}

In the example \ref Convex_hull_2/ch_graham_anderson.cpp, `ch_graham_andrew_scan()` is used to
realize Anderson's variant \cgalCite{a-readc-78} of the Graham Scan
\cgalCite{g-eadch-72}. The points are sorted counterclockwise around the leftmost
point using the `Less_rotate_ccw_2` predicate, as defined in
the concept `ConvexHullTraits_2`. According to the definition
of `Less_rotate_ccw_2`, the leftmost point is the last point in the sorted
sequence and its predecessor on the convex hull is the first point in the
sorted sequence. It is not hard to see that the preconditions of
`ch_graham_andrew_scan()` are satisfied. Anderson's variant of the
Graham scan is usually inferior to Andrew's variant because of its higher
arithmetic demand.

*/
template <class BidirectionalIterator, class OutputIterator,
class Traits>
OutputIterator
ch_graham_andrew_scan( BidirectionalIterator first,
BidirectionalIterator beyond,
OutputIterator result,
const Traits& ch_traits = Default_traits);

} /* namespace CGAL */
