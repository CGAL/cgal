
namespace CGAL {

/*!
\ingroup PkgOptimalDistances

The class `Width_default_traits_3` is a traits class for `Width_3<Traits>` 
using the three-dimensional \cgal kernel. 

\tparam K must be a model for `Kernel`.

\cgalModels `WidthTraits_3`

\sa `CGAL::Width_3<Traits>` 
\sa `WidthTraits_3` 

*/
template< typename K >
class Width_default_traits_3 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef typename K::Point_3 Point_3; 

/*!

*/ 
typedef typename K::Plane_3 Plane_3; 

/*!

*/ 
typedef typename K::Vector_3 Vector_3; 

/*!

*/ 
typedef typename K::RT RT; 

/*!

*/ 
typedef Convex_hull_traits_3<K> ChullTraits; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
Width_default_traits_3( ); 

/// @} 

/// \name Operations 
/// @{

/*!
returns the 
homogeneous \f$ x\f$-coordinate of point \f$ p\f$. 
*/ 
RT get_hx(const Point_3& p) const; 

/*!
returns the 
homogeneous \f$ y\f$-coordinate of point \f$ p\f$. 
*/ 
RT get_hy(const Point_3& p) const; 

/*!
returns the 
homogeneous \f$ z\f$-coordinate of point \f$ p\f$. 
*/ 
RT get_hz(const Point_3& p) const; 

/*!
returns the 
homogenizing coordinate of point \f$ p\f$. 
*/ 
RT get_hw(const Point_3& p) const; 

/*!
returns all homogeneous coordinates 
of point \f$ p\f$ at once. 
*/ 
void get_point_coordinates(const Point_3& p, RT& px, 
RT& py, RT& pz, RT& ph) const; 

/*!
returns the first 
coefficient of plane \f$ f\f$. 
*/ 
RT get_a(const Plane_3& f) const; 

/*!
returns the second 
coefficient of plane \f$ f\f$. 
*/ 
RT get_b(const Plane_3& f) const; 

/*!
returns the third 
coefficient of plane \f$ f\f$. 
*/ 
RT get_c(const Plane_3& f) const; 

/*!
returns the fourth 
coefficient of plane \f$ f\f$. 
*/ 
RT get_d(const Plane_3& f) const; 

/*!
returns all four plane coefficients of \f$ f\f$ at once. 
*/ 
void get_plane_coefficients(const Plane_3& f, 
RT& a, RT& b, RT& c, RT& d) 
const; 

/*!
returns a point of type 
`Point_3` with homogeneous coordinates \f$ hx\f$, \f$ hy\f$, \f$ hz\f$ and \f$ hw\f$. 
*/ 
Point_3 make_point(const RT& hx, const RT& hy, const RT& hz, 
const RT& hw) const; 

/*!
returns a plane of type `Plane_3` 
whose coefficients are \f$ a\f$, \f$ b\f$, \f$ c\f$ and \f$ d\f$. 
*/ 
Plane_3 make_plane(const RT& a, const RT& b, const 
RT& c, const RT& d) const; 

/*!
returns a vector of type `Vector_3` with the four 
coefficients \f$ a\f$, \f$ b\f$, \f$ c\f$ and 1. 
*/ 
Vector make_vector(const RT& a, const RT& b, const 
RT& c) const; 

/// @}

}; /* end Width_default_traits_3 */
} /* end namespace CGAL */
