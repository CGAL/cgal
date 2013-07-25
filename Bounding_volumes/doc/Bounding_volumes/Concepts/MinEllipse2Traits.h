
/*!
\ingroup PkgBoundingVolumesConcepts
\cgalConcept

This concept defines the requirements for traits classes of 
`CGAL::Min_ellipse_2<Traits>`. 

\cgalHasModel `CGAL::Min_ellipse_2_traits_2<K>`

\sa `CGAL::Min_ellipse_2<Traits>` 

*/

class MinEllipse2Traits {
public:

/// \name Types 
/// @{

/*!

The point type must provide default and copy constructor, 
assignment and equality test. 
*/ 
typedef unspecified_type Point; 

/*!

The ellipse type must fulfill the requirements listed below 
in the next section. 
*/ 
typedef unspecified_type Ellipse; 

/// @} 

/// \name Variables 
/// @{

/*!

The current ellipse. This variable is maintained by the algorithm, 
the user should neither access nor modify it directly. 
*/ 
Ellipse ellipse; 

/// @} 

/// \name Creation 
/// Only default and copy constructor are required.
/// @{

/*!

*/ 
Traits( ); 

/*!

*/ 
Traits( const Traits&); 

/// @}

}; /* end MinEllipse2Traits */

