
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

*/

class CircularKernel::ConstructCircularArc_2 {
public:

/// \name Operations
/// A model of this concept must provide:
/// @{

/*!
Constructs an arc from a full circle. 
*/ 
CircularKernel::Circular_arc_2 operator() 
(const CircularKernel::Circle_2 &c); 

/*!
Construct the circular arc supported by `c`, whose source is 
`p1` and whose target is `p2` when traversing the circle in 
counterclockwise direction. 
\pre `p1` and `p2` lie on `c`. 
*/ 
CircularKernel::Circular_arc_2 operator() 
(const CircularKernel::Circle_2 &c, 
const CircularKernel::Circular_arc_point_2 &p1, 
const CircularKernel::Circular_arc_point_2 &p2); 

/*!
Constructs the unique circular arc whose supporting circle is 
`c`, and whose source is the intersection of `c` 
and `c1` with index `b1`, and whose target is the intersection 
of `c` and `c2` of index `b2`, where intersections are 
ordered lexicographically, and when traversing the circle in 
counterclockwise direction. 
\pre `c` intersects both `c1` and `c2`, and the arc defined by the intersections has non-zero length. 
*/ 
CircularKernel::Circular_arc_2 operator() 
(const CircularKernel::Circle_2 &c, 
const CircularKernel::Circle_2 &c1, bool b1, 
const CircularKernel::Circle_2 &c2, bool b2); 

/*!
Same, for intersections defined by lines instead of circles. 
*/ 
CircularKernel::Circular_arc_2 operator() 
(const CircularKernel::Circle_2 &c, 
const CircularKernel::Line_2 &l1, bool b1, 
const CircularKernel::Line_2 &l2, bool b2); 

/*!
Constructs an arc that is supported by the circle of type 
`CircularKernel::Circle_2` passing through the points `p`, 
`q` and `r`. The source and target are respectively `p` 
and `r`, when traversing the supporting circle in the 
counterclockwise direction. 
Note that, depending on the orientation of the point triple 
`(p,q,r)`, `q` may not lie on the arc. 
\pre `p`, `q`, and `r` are not collinear. 
*/ 
CircularKernel::Circular_arc_2 operator() 
(const CircularKernel::Point_2 &p, 
const CircularKernel::Point_2 &q, 
const CircularKernel::Point_2 &r); 

/// @}

}; /* end CircularKernel::ConstructCircularArc_2 */

