
namespace CGAL {

/*!
\ingroup PkgCircularKernel2

\models ::CircularKernel::CircularArcPoint_2 

\sa `CGAL::Circular_arc_2<CircularKernel>`
\sa `CGAL::Line_arc_2<CircularKernel>`

*/
template< typename CircularKernel >
class Circular_arc_point_2 {
public:

/// \name Creation 
/// @{

/*! 

*/ 
Circular_arc_point_2(const CircularKernel::Point_2 &q); 

/*! 

*/ 
Circular_arc_point_2(const CircularKernel::Root_for_circles_2_2 &r); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
\f$ x\f$-coordinate of the point. 
*/ 
const CircularKernel::Root_of_2 & x(); 

/*! 
\f$ y\f$-coordinate of the point. 
*/ 
const CircularKernel::Root_of_2 & y(); 

/*! 
Returns a bounding box around the point. 
*/ 
Bbox_2 bbox() const; 

/// @}

}; /* end Circular_arc_point_2 */


/*! 
Test for equality. Two points are equal, iff their \f$ x\f$ and \f$ y\f$ coordinates are equal. 
\relates Circular_arc_point_2 
*/ 
bool operator==(const Circular_arc_point_2<CircularKernel> &p, const Circular_arc_point_2<CircularKernel> &q); 

/*! 
Test for nonequality. 
\relates Circular_arc_point_2 
*/ 
bool operator!=(const Circular_arc_point_2<CircularKernel> &p, const Circular_arc_point_2<CircularKernel> &q); 

/*! 
Returns true iff \f$ p\f$ is lexicographically smaller than \f$ q\f$, i.e. either if \f$ p.x() < q.x()\f$ 
or if \f$ p.x() == q.x()\f$ and \f$ p.y() < q.y()\f$. 
\relates Circular_arc_point_2 
*/ 
bool operator<(const Circular_arc_point_2<CircularKernel> &p, const Circular_arc_point_2<CircularKernel> &q); 

/*! 
Returns true iff \f$ p\f$ is lexicographically greater than \f$ q\f$. \relates Circular_arc_point_2 
*/ 
bool operator>(const Circular_arc_point_2<CircularKernel> &p, 
const Circular_arc_point_2<CircularKernel> &q); 

/*! 
Returns true iff \f$ p\f$ is lexicographically smaller than or equal to \f$ q\f$. 
\relates Circular_arc_point_2 
*/ 
bool operator<=(const Circular_arc_point_2<CircularKernel> &p, 
const Circular_arc_point_2<CircularKernel> &q); 

/*! 
Returns true iff \f$ p\f$ is lexicographically greater than or equal to \f$ q\f$. 
\relates Circular_arc_point_2 
*/ 
bool operator>=(const Circular_arc_point_2<CircularKernel> &p, 
const Circular_arc_point_2<CircularKernel> &q); 

/*! 

\relates Circular_arc_point_2 
*/ 
istream& operator>> (std::istream& is, Circular_arc_point_2 & cp); 

/*! 

\relates Circular_arc_point_2 
*/ 
ostream& operator<< (std::ostream& os, const Circular_arc_point_2 &ce); 

} /* end namespace CGAL */
