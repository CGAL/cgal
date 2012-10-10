
namespace CGAL {

/*!
\ingroup PkgSkinSurface3

The `Skin_surface_3` is the main class in the 
`Skin_surface_3` package. It is a model of 
the concept `SkinSurface_3`. 

The template argument must be a model of the concept 
`SkinSurfaceTraits_3`, which means that it provides predicates 
to construct a regular triangulation of the weighted points and for 
point location in the mixed complex. 

\models ::SkinSurface_3 

*/
template< typename SkinSurfaceTraits_3 >
class Skin_surface_3 {
public:

/// @}

}; /* end Skin_surface_3 */
} /* end namespace CGAL */
