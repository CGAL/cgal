namespace CGAL {

/*!
\ingroup PkgConvexHull2Functions

The function `convex_hull_2` generates the counterclockwise sequence of extreme 
points from a given set of input points. 

The default traits class `Default_traits` is the kernel in which 
the type `InputIterator::value_type` is defined. 

Requirements 
-------------- 

<OL> 
<LI>`InputIterator::value_type` and `OutputIterator::value_type` 
are equivalent to `Traits::Point_2`. 
<LI>`Traits` contains the following subset of types from 
the concept `ConvexHullTraits_2` and their corresponding member 
functions that return instances of these types: 
<UL> 
<LI>`Traits::Point_2`, 
<LI>`Traits::Less_signed_distance_to_line_2`, 
<LI>`Traits::Equal_2`, 
<LI>`Traits::Less_xy_2`, 
<LI>`Traits::Less_yx_2`, 
<LI>`Traits::Left_turn_2`. 
</UL> 
</OL> 

\sa `CGAL::ch_akl_toussaint` 
\sa `CGAL::ch_bykat` 
\sa `CGAL::ch_eddy` 
\sa `CGAL::ch_graham_andrew` 
\sa `CGAL::ch_jarvis` 
\sa `CGAL::ch_melkman` 

Implementation 
-------------- 

One of two algorithms is used, 
depending on the type of iterator used to specify the input points. For 
input iterators, the algorithm used is that of Bykat \cite b-chfsp-78, which 
has a worst-case running time of \f$ O(n h)\f$, where \f$ n\f$ is the number of input 
points and \f$ h\f$ is the number of extreme points. For all other types of 
iterators, the \f$ O(n \log n)\f$ algorithm of of Akl and Toussaint 
\cite at-fcha-78 is used. 

generates the counterclockwise sequence of extreme points
of the points in the range [`first`,`beyond`).
The resulting sequence is placed starting at position
`result`, and the past-the-end iterator for the resulting
sequence is returned. It is not specified at which point the
cyclic sequence of extreme points is cut into a linear sequence.
\pre The source range [`first`,`beyond`) does not contain `result`.
*/
template <class InputIterator, class OutputIterator>
OutputIterator
convex_hull_2(InputIterator first, InputIterator beyond,
OutputIterator result,
const Traits & ch_traits = Default_traits );

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgConvexHull2Subsequence

The function `lower_hull_points_2` generates the counterclockwise sequence of extreme 
points on the lower hull of a given set of input points. 

The default traits class `Default_traits` is the kernel in which the 
type `InputIterator::value_type` is defined. 

The different treatment by `CGAL::upper_hull_points_2` of the case that 
all points are equal ensures that concatenation of lower and upper hull 
points gives the sequence of extreme points. 

Requirements 
-------------- 

<OL> 
<LI>`InputIterator::value_type` and `OutputIterator::value_type` 
are equivalent to `Traits::Point_2`. 
<LI>`Traits` contains the following subset of types from 
the concept `ConvexHullTraits_2` and their corresponding member 
functions that return instances of these types: 
<UL> 
<LI>`Traits::Point_2`, 
<LI>`Traits::Equal_2`, 
<LI>`Traits::Less_xy_2`, 
<LI>`Traits::Less_yx_2`, 
<LI>`Traits::Left_turn_2`. 
</UL> 
</OL> 

\sa `CGAL::ch_graham_andrew` 
\sa `CGAL::ch_graham_andrew_scan` 
\sa `CGAL::upper_hull_points_2` 

Implementation 
-------------- 

This function uses Andrew's variant of Graham's scan algorithm 
\cite a-aeach-79, \cite m-mdscg-84. The algorithm has worst-case running time 
of \f$ O(n \log n)\f$ for \f$ n\f$ input points. 

generates the counterclockwise sequence of extreme points
on the lower hull of the points in the range [`first`,
`beyond`). The resulting sequence is placed starting at
position `result`, and the past-the-end iterator for
the resulting sequence is returned.
The sequence starts with the leftmost point;
the rightmost point is not included.
If there is only one extreme point (<I>i.e.</I>, leftmost and
rightmost point are equal) the extreme point is reported.
\pre The source range [`first`,`beyond`) does not contain
`result`.
*/
template <class InputIterator, class OutputIterator>
OutputIterator
lower_hull_points_2(InputIterator first, InputIterator beyond,
OutputIterator result,
const Traits & ch_traits = Default_traits );

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgConvexHull2Subsequence

The function `upper_hull_points_2` generates the counterclockwise sequence of extreme 
points on the upper hull of a given set of input points. 

The default traits class `Default_traits` is the kernel in which the 
type `InputIterator::value_type` is defined. 

The different treatment by `CGAL::lower_hull_points_2` of the case that 
all points are equal ensures that concatenation of lower and upper hull 
points gives the sequence of extreme points. 

Requirements 
-------------- 

<OL> 
<LI>`InputIterator::value_type` and `OutputIterator::value_type` 
are equivalent to `Traits::Point_2`. 
<LI>`Traits` contains the following subset of types from 
the concept `ConvexHullTraits_2` and their corresponding member 
functions that return instances of these types: 
<UL> 
<LI>`Traits::Point_2`, 
<LI>`Traits::Equal_2`, 
<LI>`Traits::Less_xy_2`, 
<LI>`Traits::Less_yx_2`, 
<LI>`Traits::Left_turn_2`. 
</UL> 
</OL> 

\sa `CGAL::ch_graham_andrew` 
\sa `CGAL::ch_graham_andrew_scan` 
\sa `CGAL::lower_hull_points_2` 

Implementation 
-------------- 

This function uses Andrew's 
variant of Graham's scan algorithm \cite a-aeach-79, \cite m-mdscg-84. The algorithm 
has worst-case running time of \f$ O(n \log n)\f$ for \f$ n\f$ input points. 

generates the counterclockwise sequence of extreme points
on the upper hull of the points in the range [`first`,
`beyond`). The resulting sequence is placed starting at
position `result`, and the past-the-end iterator for
the resulting sequence is returned.
The sequence starts with the rightmost point,
the leftmost point is not included.
If there is only one extreme point (<I>i.e.</I>, the leftmost and
rightmost point are equal), the extreme point is not reported.
\pre The source range [`first`,`beyond`) does not contain
`result`.
*/
template <class InputIterator, class OutputIterator>
OutputIterator
upper_hull_points_2(InputIterator first, InputIterator beyond,
OutputIterator result,
const Traits & ch_traits = Default_traits );

} /* namespace CGAL */
