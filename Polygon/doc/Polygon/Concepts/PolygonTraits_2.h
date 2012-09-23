
/*!
\ingroup PkgPolygon2Concepts
\cgalconcept

The `CGAL::Polygon_2` class and the functions that implement the
functionality found in that class each are parameterized by a traits
class that defines the primitives used in the algorithms.  The concept
`PolygonTraits_2` defines this common set of requirements.

The requirements of `PolygonTraits_2` are a subset of the kernel
requirements.  We only list the types and methods which are required
and refer to the description of the kernel concept for details.

\hasModel The kernels supplied by \cgal are models of `PolygonTraits_2`. 
\hasModel `CGAL::Projection_traits_xy_3<K>`
\hasModel `CGAL::Projection_traits_yz_3<K>`
\hasModel `CGAL::Projection_traits_zx_3<K>`

\sa `CGAL::Polygon_2<PolygonTraits_2, Container>`

*/

class PolygonTraits_2 {
public:

/// \name Creation
/// @{
  /// 
  PolygonTraits_2();
  /// 
  PolygonTraits_2(const PolygonTraits_2&);
/// @}

/// \name Types 
/// @{

/*! 

*/ 
typedef Hidden_type FT; 

/*! 
The point type. 
*/ 
typedef Hidden_type Point_2; 

/*! 
The segment type. 
*/ 
typedef Hidden_type Segment_2; 

/*! 

*/ 
typedef Hidden_type Construct_segment_2; 

/*! 

*/ 
typedef Hidden_type Equal_2; 

/*! 

*/ 
typedef Hidden_type Less_xy_2; 

/*! 

*/ 
typedef Hidden_type Less_yx_2; 

/*! 

*/ 
typedef Hidden_type Compare_x_2; 

/*! 

*/ 
typedef Hidden_type Compare_y_2; 

/*! 

*/ 
typedef Hidden_type Orientation_2; 

/*! 
Computes the signed area of the oriented 
triangle defined by 3 `Point_2` passed as arguments. 
*/ 
typedef Hidden_type Compute_area_2; 

/// @} 

/// \name Operations 
/// The following functions that create instances of the above predicate object types must exist.
/// @{

/*! 

*/ 
Equal_2 equal_2_object(); 

/*! 

*/ 
Less_xy_2 less_xy_2_object(); 

/*! 

*/ 
Less_yx_2 less_yx_2_object(); 

/*! 

*/ 
Compare_y_2 compare_y_2_object(); 

/*! 

*/ 
Compare_x_2 compare_x_2_object(); 

/*! 

*/ 
Orientation_2 orientation_2_object(); 

/*! 

*/ 
Compute_area_2 compute_area_2_object(); 

/*! 

*/ 
Construct_segment_2 construct_segment_2_object(); 

/// @}

}; /* end PolygonTraits_2 */

