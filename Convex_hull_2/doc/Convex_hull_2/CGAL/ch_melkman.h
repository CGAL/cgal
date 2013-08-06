namespace CGAL {

/*!
\ingroup PkgConvexHull2Functions

generates the counterclockwise sequence of extreme points
of the points in the range [`first`, `beyond`). 
The resulting sequence is placed starting at
position `result`, and the past-the-end iterator for
the resulting sequence is returned.
\pre The source range [`first`,`beyond`) corresponds to a simple polyline. [`first`,`beyond`) does not contain `result`.

The default traits class `Default_traits` is the kernel in which the 
value type of `InputIterator` is defined. 

\cgalHeading{Requirements}

<OL> 
<LI>The value type of `InputIterator` and `OutputIterator` 
is equivalent to `Traits::Point_2`. 
<LI>`Traits` contains the following subset of types from 
the concept `ConvexHullTraits_2` and their corresponding member 
functions that return instances of these types: 
<UL> 
<LI>`Traits::Point_2`, 
<LI>`Traits::Equal_2`, 
<LI>`Traits::Less_xy_2`, 
<LI>`Traits::Left_turn_2`. 
</UL> 
</OL> 

\sa `CGAL::ch_akl_toussaint()` 
\sa `CGAL::ch_bykat()` 
\sa `CGAL::ch_eddy()` 
\sa `CGAL::ch_graham_andrew()` 
\sa `CGAL::ch_jarvis()` 
\sa `CGAL::ch_melkman()` 
\sa `CGAL::convex_hull_2()` 

\cgalHeading{Implementation}

It uses an implementation of Melkman's algorithm \cgalCite{m-olcch-87}. 
Running time of this is linear.


*/
template <class InputIterator, class OutputIterator>
OutputIterator
ch_melkman( InputIterator first, InputIterator last, 
OutputIterator result, 
const Traits& ch_traits = Default_traits);

} /* namespace CGAL */
