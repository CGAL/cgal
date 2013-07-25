
namespace CGAL {

/*!
\ingroup PkgSurfaceMesher3Classes

A 3D gray image is a 
tri-dimensional array that associates a scalar value to each triple of 
integer \f$ (x, y, z)\f$ in the range of the image. A trilinear interpolation 
algorithm provides a map \f$ f : \mathbb{R}^3 \longrightarrow \mathbb{R}\f$. 

The class `Gray_level_image_3` is a 3D gray image loader and a model 
of the concept `ImplicitFunction`. 
An object of the class `Gray_level_image_3` is created with a parameter 
`iso_value` and then its `operator()` implements 
the function `sign of (f(p) - iso)`, for \f$ p \in \mathbb{R}^3\f$. 
Plugging such a function in the creation of the `Implicit_surface_3` 
object given as parameter to `make_surface_mesh()` yields 
a mesh approximating the level with value `iso` 
in the input 3D gray image. 

`Gray_level_image_3` provides an interface with an auxiliary library called 
<I>CGAL_ImageIO</I>. An executable that uses `Gray_level_image_3` must be linked with 
the <I>CGAL_ImageIO</I> library.

The library <I>CGAL_ImageIO</I> and therefore `Gray_level_image_3` support 
several types of 3D images: INRIMAGE (extension .inr[.gz]), GIS (extension 
.dim, of .ima[.gz]), and ANALYZE (extension .hdr, or .img[.gz]). 

\cgalModels `ImplicitFunction`

\sa `ImplicitFunction`
\sa `Implicit_surface_3<Traits, Function>`
\sa `make_surface_mesh` 

*/
template< typename FT_, typename Point_ >
class Gray_level_image_3 {
public:

/// \name Types 
/// @{

/*!
the numerical type `FT`.
*/ 
typedef FT_ FT; 

/*!
the point type. 
*/ 
typedef Point_ Point; 


/// @} 

/// \name Creation 
/// @{

/*!
`filename` is the path to a file of a type supported by <I>ImageIO</I>. 

`iso_value` is an isovalue of the interpolation function \f$ f\f$. 

*/ 
Gray_level_image_3(const char* filename, FT iso_value); 

/// @}


/// \name Operations
/// @{

/*! Returns the sign of  \f$ f(p)\f$  - `iso_value`. 
 */
  FT operator(const Point& p) const;
/// @
}; /* end Gray_level_image_3 */
} /* end namespace CGAL */
