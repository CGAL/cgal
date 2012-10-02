
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3

The class `Fixed_alpha_shape_vertex_base_3` is the default model for the concept 
`FixedAlphaShapeVertex_3`. 

The class has two parameters : the traits class `Traits` 
which provides the type for the points or the weighted points. 
The second parameter `Vb` is a base class instantiated by default 
with `CGAL::Triangulation_vertex_base_3<Traits>`. 

\models ::FixedAlphaShapeVertex_3 

*/
template< typename Traits, typename Vb >
class Fixed_alpha_shape_vertex_base_3 : public Vb {
public:

/// @}

}; /* end Fixed_alpha_shape_vertex_base_3 */
} /* end namespace CGAL */
