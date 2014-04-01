
/*!
\ingroup PkgSurfaceMesher3Concepts
\cgalConcept

The concept `ImplicitFunction` describes a function object 
whose `operator()` computes the values of a function 
\f$ f : \mathbb{R}^3 \longrightarrow \mathbb{R}\f$. 

\cgalHasModel CGAL::Gray_level_image_3
\cgalHasModel any pointer to a function of type `FT (*)(Point)`

\sa `CGAL::Implicit_surface_3<Traits, Function>`
\sa `CGAL::make_surface_mesh()` 

*/

class ImplicitFunction {
public:

/// \name Types 
/// @{
///The following types aren't required for any pointer to a function of type `FT (*)(Point)`.
/*!
Number type 
*/ 
typedef unspecified_type FT; 

/*!
Point type 
*/ 
typedef unspecified_type Point; 

/// @} 

/// \name Operations 
/// @{

/*!
Returns the value \f$ f(p)\f$, where \f$ p \in\mathbb{R}^3\f$. 
*/ 
FT operator()(Point p); 

/// @}

}; /* end ImplicitFunction */

