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
value type of  `ForwardIterator` is defined.

\cgalHeading{Requirements}

<OL>
<LI>The value type of `InputIterator` and
`OutputIterator` is equivalent to `Traits::Point_2`.
<LI>`Traits` defines the following subset of types from
the concept `ConvexHullTraits_2` and their corresponding member
functions that return instances of these types:
<UL>
<LI>`Traits::Point_2`,
<LI>`Traits::Less_signed_distance_to_line_2`,
<LI>`Traits::Left_turn_2`,
<LI>`Traits::Less_xy_2`,
<LI>`Traits::Equal_2`.
</UL>
</OL>

\sa `CGAL::ch_akl_toussaint()`
\sa `CGAL::ch_eddy()`
\sa `CGAL::ch_graham_andrew()`
\sa `CGAL::ch_jarvis()`
\sa `CGAL::ch_melkman()`
\sa `CGAL::convex_hull_2()`

\cgalHeading{Implementation}

This function implements the non-recursive variation of
Eddy's algorithm \cgalCite{e-nchap-77} described in \cgalCite{b-chfsp-78}.
This algorithm requires \f$ O(n h)\f$ time
in the worst case for \f$ n\f$ input points with \f$ h\f$ extreme points.


*/
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_bykat( InputIterator first,
InputIterator beyond,
OutputIterator result,
const Traits & ch_traits = Default_traits);

} /* namespace CGAL */
