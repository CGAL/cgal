
/*!
\ingroup PkgSurfaceMesher3Concepts
\cgalconcept

The concept `ImplicitFunction` describes a function object 
whose `operator()` computes the values of a function 
\f$ f : \R^3 \longrightarrow \R\f$. 

\hasModel `CGAL::Gray_level_image_3`
\hasModel any pointer to a function of type `FT (*)(Point)`

\sa `CGAL::Implicit_surface_3<Traits, Function>`
\sa `CGAL::make_surface_mesh` 

*/

class ImplicitFunction {
public:

/// \name Types 
/// @{

/*! 
Number type 
*/ 
typedef Hidden_type FT; 

/*! 
Point type 
*/ 
typedef Hidden_type Point; 

/// @} 

/// \name Operations 
/// @{

/*! 
Returns the value \f$ f(p)\f$, where \f$ p \in\R^3\f$. 
*/ 
FT operator()(Point p); 

/// @}

}; /* end ImplicitFunction */

