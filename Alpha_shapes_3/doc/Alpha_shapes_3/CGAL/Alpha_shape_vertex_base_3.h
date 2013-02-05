
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3

The class `Alpha_shape_vertex_base_3` is the default model for the concept 
`AlphaShapeVertex_3`. 

\tparam Traits must be a geometric traits class
which provides the type for the points or the weighted points. 
\tparam Vb must be a vertex base class instantiated by default 
with `CGAL::Triangulation_vertex_base_3<Traits>`. 
\tparam ExactAlphaComparisonTag is a tag that, when set to 
`Tag_true`, triggers exact comparisons between alpha values. See the description 
provided in the documentation of `Alpha_shape_3` for more details. The default value is `Tag_false`. 
\tparam WeightedTag is used only if `ExactAlphaComparisonTag` is `Tag_true`. It 
must be `Tag_true` if the underlying triangulation of the alpha shape to be used is a regular triangulation 
and `Tag_false` otherwise. The default is `Tag_false`. 

\cgalModels `AlphaShapeVertex_3`

*/
template< typename Traits, typename Vb, typename ExactAlphaComparisonTag, typename WeightedTag >
class Alpha_shape_vertex_base_3 : public Vb {
public:

/// @}

}; /* end Alpha_shape_vertex_base_3 */
} /* end namespace CGAL */
