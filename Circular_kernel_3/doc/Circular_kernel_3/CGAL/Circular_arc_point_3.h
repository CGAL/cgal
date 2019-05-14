
namespace CGAL {

/*!
\ingroup PkgSphericalKernel3GeometricClasses

\cgalModels `SphericalKernel::CircularArcPoint_3`

\sa `CGAL::Circular_arc_3<SphericalKernel>`
\sa `CGAL::Line_arc_3<SphericalKernel>`

*/
template< typename SphericalKernel >
class Circular_arc_point_3 {
public:

/// \name Creation 
/// @{

/*!

*/ 
Circular_arc_point_3(const Point_3<SphericalKernel> &q); 

/*!

*/ 
Circular_arc_point_3(const SphericalKernel::Root_for_spheres_2_3 &r); 

/// @} 

/// \name Access Functions 
/// @{

/*!
\f$ x\f$-coordinate of the point. 
*/ 
const SphericalKernel::Root_of_2 & x(); 

/*!
\f$ y\f$-coordinate of the point. 
*/ 
const SphericalKernel::Root_of_2 & y(); 

/*!
\f$ z\f$-coordinate of the point. 
*/ 
const SphericalKernel::Root_of_2 & z(); 

/*!
Returns a bounding box around the point. 
*/ 
Bbox_3 bbox() const; 

/// @}

}; /* end Circular_arc_point_3 */


/*!
Test for equality. Two points are equal, iff their \f$ x\f$, \f$ y\f$ and \f$ z\f$ 
coordinates are equal. 
\relates Circular_arc_point_3 
*/ 
bool operator==(const Circular_arc_point_3<SphericalKernel> &p, const Circular_arc_point_3<SphericalKernel> &q); 

/*!
Test for nonequality. 
\relates Circular_arc_point_3 
*/ 
bool operator!=(const Circular_arc_point_3<SphericalKernel> &p, const Circular_arc_point_3<SphericalKernel> &q); 

/*!
Returns true iff `p` is lexicographically smaller than `q`, i.e.\ either 
if `p.x() < q.x()` or if `p.x() == q.x()` and `p.y() < q.y()`
or if `p.x() == q.x()` and `p.y() == q.y()` and `p.z() < q.z()`. 
\relates Circular_arc_point_3 
*/ 
bool operator<(const Circular_arc_point_3<SphericalKernel> &p, const Circular_arc_point_3<SphericalKernel> &q); 

/*!
Returns true iff `p` is lexicographically greater than `q`. 
\relates Circular_arc_point_3 
*/ 
bool operator>(const Circular_arc_point_3<SphericalKernel> &p, const Circular_arc_point_3<SphericalKernel> &q); 

/*!
Returns true iff `p` is lexicographically smaller than or equal to `q`. 
\relates Circular_arc_point_3 
*/ 
bool operator<=(const Circular_arc_point_3<SphericalKernel> &p, const Circular_arc_point_3<SphericalKernel> &q); 

/*!
Returns true iff `p` is lexicographically greater than or equal to `q`. 
\relates Circular_arc_point_3 
*/ 
bool operator>=(const Circular_arc_point_3<SphericalKernel> &p, const Circular_arc_point_3<SphericalKernel> &q); 

/*!

\relates Circular_arc_point_3 
*/ 
istream& operator>> (std::istream& is, Circular_arc_point_3 & p); 

/*!

\relates Circular_arc_point_3 
*/ 
ostream& operator<< (std::ostream& os, const Circular_arc_point_3 &p); 

} /* end namespace CGAL */
