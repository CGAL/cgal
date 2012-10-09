namespace CGAL {

/*!
\ingroup PkgConvexHull2Traits

\deprecated The functionality of this class has been generalized to other packages than 2D convex hulls. 
The more general class `Projection_traits_yz_3` can be found in the 2D and 3D Linear Geometric Kernel. 
Note that the deprecated class was templated by a point class, whereas the new class 
is templated by a geometric kernel. 

The class `Convex_hull_projective_yz_traits_2` serves as a traits class for all the two-dimensional 
convex hull and extreme point calculation function. This class can be 
used to compute the convex hull of a set of 3D points projected onto the 
\f$ yz\f$ plane (<I>i.e.</I>, by ignoring the \f$ x\f$ coordinate). 

\models ::ConvexHullTraits_2 

\sa `CGAL::Convex_hull_constructive_traits_2<R>` 
\sa `CGAL::Convex_hull_projective_xz_traits_2<Point_3>` 
\sa `CGAL::Convex_hull_projective_xy_traits_2<Point_3>` 
\sa `CGAL::Convex_hull_traits_2<R>` 

*/
template< typename Point_3 >
class Convex_hull_projective_yz_traits_2 {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef Point_3 Point_2; 

/*! 

*/ 
typedef Less_xy_plane_yz_2<Point_3> Less_xy_2; 

/*! 

*/ 
typedef Less_yx_plane_yz_2<Point_3> Less_yx_2; 

/*! 

*/ 
typedef Less_dist_to_line_plane_yz_2<Point_3> 
Less_signed_distance_to_line_2; 

/*! 

*/ 
typedef Less_rotate_ccw_plane_yz_2<Point_3> Less_rotate_ccw_2; 

/*! 

*/ 
typedef Left_turn_plane_yz_2<Point_3> Left_turn_2; 

/*! 

*/ 
typedef Equal_xy_plane_yz_2<Point_3> Equal_2; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Convex_hull_projective_yz_traits_2(); 

/// @} 

/// \name Operations 
/// @{

/*! 

*/ 
Less_xy_2 less_xy_2_object(); 

/*! 

*/ 
Less_yx_2 less_yx_2_object(); 

/*! 

*/ 
Less_signed_distance_to_line_2 
less_signed_distance_to_line_2_object(); 

/*! 

*/ 
Less_rotate_ccw_2 
less_rotate_ccw_2_object(); 

/*! 

*/ 
Left_turn_2 left_turn_2_object(); 

/*! 

*/ 
Equal_2 equal_2_object(); 

/// @}

}; /* end Convex_hull_projective_yz_traits_2 */
} /* end namespace CGAL */
