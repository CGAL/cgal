namespace CGAL {

/*!

The class `Random_convex_set_traits_2` serves as a traits class 
for the function `random_convex_set_2()`. 

\cgalModels `RandomConvexSetTraits_2`

*/
template< typename Kernel >
class Random_convex_set_traits_2 {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef Kernel::Point_2 Point_2; 

/*!

*/ 
typedef Kernel::FT FT; 

/*!
function object class derived from 
`std::binary_function<Point_2, Point_2, Point_2>` 

*/ 
typedef unspecified_type Sum; 

/*!
function object class derived from 
`std::binary_function<Point_2, Point_2, Point_2>` 
*/ 
typedef unspecified_type Scale; 

/*!
function object class derived from 
`std::unary_function<Point_2, FT>` 
*/ 
typedef unspecified_type Max_coordinate; 

/*!
function object class derived from 
`std::binary_function<Point_2, Point_2, bool>` 
*/ 
typedef unspecified_type Angle_less; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor 
*/ 
Random_convex_set_traits_2(); 

/// @} 

/// \name Operations 
/// @{

/*!
returns CGAL::ORIGIN. 
*/ 
Point_2 origin() const; 

/// @}

}; /* end Random_convex_set_traits_2 */
} /* end namespace CGAL */
