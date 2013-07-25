
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

*/

class CircularKernel::ConstructLineArc_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Constructs the line segment supported by `l`, whose source 
is `p1` and whose target is `p2`. 
\pre `p1` and `p2` lie on `l`. 
*/ 
CircularKernel::Line_arc_2 operator() 
(const CircularKernel::Line_2 &l, 
const CircularKernel::Circular_arc_point_2 &p1, 
const CircularKernel::Circular_arc_point_2 &p2); 

/*!

*/ 
CircularKernel::Line_arc_2 operator() 
(const CircularKernel::Segment_2 &s); 

/*!

*/ 
CircularKernel::Line_arc_2 operator() 
(const CircularKernel::Point_2 &p1, 
const CircularKernel::Point_2 &p2); 

/*!
Constructs the line segment whose supporting line is `l`, whose 
source endpoint is the \f$ b_1^{th}\f$ intersection of `l` with `c1`, 
and whose target endpoint is the \f$ b_2^{th}\f$ intersection of `l` 
and `c2`, where intersections are ordered lexicographically. 
\pre `l` intersects both `c1` and `c2`, and the arc defined by the intersections has non-zero length. 
*/ 
CircularKernel::Line_arc_2 operator() 
(const CircularKernel::Line_2 &l, 
const CircularKernel::Circle_2 &c1, bool b1, 
const CircularKernel::Circle_2 &c2, bool b2); 

/*!
Same, for intersections defined by lines instead of circles. 
*/ 
CircularKernel::Line_arc_2 operator() 
(const CircularKernel::Line_2 &l, 
const CircularKernel::Line_2 &l1, 
const CircularKernel::Line_2 &l2); 

/// @}

}; /* end CircularKernel::ConstructLineArc_2 */

