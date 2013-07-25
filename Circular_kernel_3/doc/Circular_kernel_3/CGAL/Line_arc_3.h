
namespace CGAL {

/*!
\ingroup PkgSphericalKernel3GeometricClasses

\cgalModels `SphericalKernel::LineArc_3`

\sa `CGAL::Circular_arc_point_3<SphericalKernel>`
\sa `CGAL::Circular_arc_3<SphericalKernel>`
*/
template< typename SphericalKernel >
class Line_arc_3 {
public:

/// \name Creation 
/// @{

/*!
Construct the line segment supported by `l`, whose source 
is `p1`, and whose target is `p2`. 
\pre `p1` and `p2` lie on `l`. `p1` and `p2` are different. 
*/ 
Line_arc_3(const Line_3<SphericalKernel> &l, 
const Circular_arc_point_3<SphericalKernel> &p1, 
const Circular_arc_point_3<SphericalKernel> &p2); 

/*!
Same. 
*/ 
Line_arc_3(const Line_3<SphericalKernel> &l, 
const Point_3<SphericalKernel> &p1, 
const Point_3<SphericalKernel> &p2); 

/*!

*/ 
Line_arc_3(const Segment_3<SphericalKernel> &s); 

/// @} 

/// \name Access Functions 
/// @{

/*!

*/ 
Line_3<SphericalKernel> supporting_line(); 

/*!

*/ 
Circular_arc_point_3<SphericalKernel> source(); 

/*!

*/ 
Circular_arc_point_3<SphericalKernel> target(); 

/*!
Constructs the minimum vertex according to the lexicographic ordering 
of coordinates. 
*/ 
Circular_arc_point_3<SphericalKernel> min(); 

/*!
Same for the maximum vertex. 
*/ 
Circular_arc_point_3<SphericalKernel> max(); 

/// @} 

/// \name Query Functions 
/// @{

/*!
Returns true `iff` the segment is 
vertical. 
*/ 
bool is_vertical(); 

/*!
Test for equality. Two segments are equal, iff their non-oriented 
supporting lines are equal (i.e.\ they define the same set of 
points), and their endpoints are the same. 
*/ 
bool operator==(const Line_arc_3<SphericalKernel> &s1, const Line_arc_3<SphericalKernel> &s2); 

/*!
Test for nonequality. 
*/ 
bool operator!=(const Line_arc_3<SphericalKernel> &s1, const Line_arc_3<SphericalKernel> &s2); 

/*!
The format for input/output is, for each line arc: a `Line_3` 
(the supporting line) and two `Circular_arc_point_3` (the two endpoints), 
under the condition that the endpoints are actually lying on the line. 
*/ 
istream& operator>> (std::istream& is, Line_arc_3 & ca); 

/*!
The format for input/output is, for each line arc: a `Line_3` 
(the supporting line) and two `Circular_arc_point_3` (the two endpoints), 
under the condition that the endpoints are actually lying on the line. 
*/ 
ostream& operator<< (std::ostream& os, const Line_arc_3 & ca); 

/// @}

}; /* end Line_arc_3 */
} /* end namespace CGAL */
