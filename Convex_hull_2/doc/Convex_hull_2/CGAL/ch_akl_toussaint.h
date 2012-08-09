namespace CGAL {

/*!
\ingroup PkgConvexHull2Functions

The function `ch_akl_toussaint` generates the counterclockwise sequence of extreme 
points from a given set of input points. 

The default traits class `Default_traits` is the kernel in which the 
type `ForwardIterator::value_type` is defined. 

Requirements 
-------------- 

<OL> 
<LI>`ForwardIterator::value_type` and 
`OutputIterator::value_type` 
are equivalent to `Traits::Point_2`. 
<LI>`Traits` defines the following subset of types from 
the concept `ConvexHullTraits_2` and their corresponding member 
functions that return instances of these types: 
<UL> 
<LI>`Traits::Point_2`, 
<LI>`Traits::Less_xy_2`, 
<LI>`Traits::Less_yx_2`, 
<LI>`Traits::Left_turn_2`, 
<LI>`Traits::Equal_2`. 
</UL> 
</OL> 

\sa `CGAL::ch_bykat` 
\sa `CGAL::ch_eddy` 
\sa `CGAL::ch_graham_andrew` 
\sa `CGAL::ch_jarvis` 
\sa `CGAL::ch_melkman` 
\sa `CGAL::convex_hull_2` 

Implementation 
-------------- 

This function uses the algorithm of Akl and 
Toussaint \cite at-fcha-78 that requires \f$ O(n \log n)\f$ time for \f$ n\f$ input 
points. 

generates the counterclockwise sequence of extreme points
of the points in the range [`first`,`beyond`).
The resulting sequence is placed starting at position
`result`, and the past-the-end iterator for the resulting
sequence is returned. It is not specified at which point the
cyclic sequence of extreme points is cut into a linear sequence.
\pre The source range [`first`,`beyond`) does not contain `result`.
*/
template <class ForwardIterator, class OutputIterator, class Traits>
OutputIterator
ch_akl_toussaint(ForwardIterator first, ForwardIterator beyond,
OutputIterator result,
const Traits& ch_traits = Default_traits());

} /* namespace CGAL */
