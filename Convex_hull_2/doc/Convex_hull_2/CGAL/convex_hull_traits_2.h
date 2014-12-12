namespace CGAL {

/*!
\ingroup PkgConvexHull2Traits

The class `Convex_hull_traits_2` serves as a traits class for all the two-dimensional 
convex hull and extreme point calculation function. This class corresponds 
to the default traits class for these functions. 

\cgalModels `ConvexHullTraits_2`

\sa `CGAL::Convex_hull_constructive_traits_2<R>` 
\sa `CGAL::Projection_traits_xy_3<K>`
\sa `CGAL::Projection_traits_yz_3<K>`
\sa `CGAL::Projection_traits_xz_3<K>`

*/
template< typename R >
class Convex_hull_traits_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef R::Point_2 Point_2; 

/*!

*/ 
typedef R::Less_xy Less_xy_2; 

/*!

*/ 
typedef R::Less_yx Less_yx_2; 

/*!

*/ 
typedef R::Less_signed_distance_to_line_2 
Less_signed_distance_to_line_2; 

/*!

*/ 
typedef R::Less_rotate_ccw_2 Less_rotate_ccw_2; 

/*!

*/ 
typedef R::Left_turn_2 Left_turn_2; 

/*!

*/ 
typedef R::Equal_2 Equal_2; 

/// @} 

/// \name Creation 
/// @{

/*!
copy constructor. 
*/ 
Convex_hull_traits_2(Convex_hull_traits_2& t); 

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
Less_rotate_ccw_2 less_rotate_ccw_2_object(); 

/*!

*/ 
Left_turn_2 left_turn_2_object(); 

/*!

*/ 
Equal_2 equal_2_object(); 

/// @}

}; /* end Convex_hull_traits_2 */
} /* end namespace CGAL */
