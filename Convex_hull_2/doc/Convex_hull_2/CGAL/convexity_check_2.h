namespace CGAL {

/*!
\ingroup PkgConvexHull2Convexity

The function `is_ccw_strongly_convex_2()` determines if a given sequence of points defines 
a counterclockwise-oriented, strongly convex polygon. 
It returns `true`, iff the point elements in 
[`first`,`beyond`)
form a counterclockwise-oriented strongly convex polygon.

A set of points is said to be strongly convex 
if it consists of only extreme points 
(<I>i.e.</I>, vertices of the convex hull). 

The default traits class `Default_traits` is the kernel in which the 
value type of `ForwardIterator` is defined. 

\cgalHeading{Requirements}

`Traits` contains the following subset of types from 
the concept `ConvexHullTraits_2` and their corresponding member 
functions that return instances of these types: 
<UL> 
<LI>`Traits::Less_xy_2`, 
<LI>`Traits::Equal_2`, 
<LI>`Traits::Left_turn_2`. 
</UL> 

\sa `CGAL::is_cw_strongly_convex_2()` 
\sa `CGAL::is_strongly_convex_3()` 

\cgalHeading{Implementation}

The algorithm requires \f$ O(n)\f$ time for a set of \f$ n\f$ input points. 



*/
template <class ForwardIterator, class Traits>
bool
is_ccw_strongly_convex_2( 
ForwardIterator first,
ForwardIterator beyond,
const Traits & ch_traits = Default_traits);

} /* namespace CGAL */

namespace CGAL {

/*!
\ingroup PkgConvexHull2Convexity

The function `is_cw_strongly_convex_2()` determines if a given sequence of points defines 
a clockwise-oriented, strongly convex polygon. 
It returns `true`, iff the point elements in 
[`first`,`beyond`)
form a clockwise-oriented strongly convex polygon.

A set of points is said to be strongly convex 
<A NAME="Index_anchor_76"></A> 
if it consists of only extreme points 
(<I>i.e.</I>, vertices of the convex hull). 

The default traits class `Default_traits` is the kernel in which the 
value type of `ForwardIterator` is defined. 

\cgalHeading{Requirements}

`Traits` contains the following subset of types from 
the concept `ConvexHullTraits_2` and their corresponding member 
functions that return instances of these types: 
<UL> 
<LI>`Traits::Equal_2`, 
<LI>`Traits::Less_xy_2`, 
<LI>`Traits::Left_turn_2`. 
</UL> 

\sa `CGAL::is_ccw_strongly_convex_2()` 
\sa `CGAL::is_strongly_convex_3()` 

\cgalHeading{Implementation}

The algorithm requires \f$ O(n)\f$ time for a set of \f$ n\f$ input points. 



*/
template <class ForwardIterator, class Traits>
bool
is_cw_strongly_convex_2( 
ForwardIterator first,
ForwardIterator beyond,
const Traits & ch_traits = Default_traits);

} /* namespace CGAL */
