
namespace CGAL {

/*!
\ingroup PkgSkinSurface3

The `Union_of_balls_3` is used to represent a skin surface with shrink 
factor equal to one, which is the boundary of the union of the input 
balls. This case is handled separately since the mixed complex is 
equal to the power diagram, which has a much simpler structure. 

The template argument must be a model of the concept 
`SkinSurfaceTraits_3`, which means that it provides the 
predicates to construct a regular triangulation of the weighted 
points. 

\cgalModels `SkinSurface_3`

*/
template< typename SkinSurfaceTraits_3 >
class Union_of_balls_3 {
public:

/// @}

}; /* end Union_of_balls_3 */
} /* end namespace CGAL */
