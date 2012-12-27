
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3

The class `Fixed_alpha_shape_vertex_base_3` is the default model for the concept 
`FixedAlphaShapeVertex_3`. 

\tparam Traits must provide the type for the points or the weighted points. 
\tparam Vb must be a vertex base class and it is instantiated by default 
with `Triangulation_vertex_base_3<Traits>`. 

\cgalModels `FixedAlphaShapeVertex_3`

*/
template< typename Traits, typename Vb >
class Fixed_alpha_shape_vertex_base_3 : public Vb {
public:

/// @}

}; /* end Fixed_alpha_shape_vertex_base_3 */
} /* end namespace CGAL */
