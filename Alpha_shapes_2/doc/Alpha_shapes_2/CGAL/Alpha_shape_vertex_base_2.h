
namespace CGAL {

/*!
\ingroup PkgAlphaShape2

The class `Alpha_shape_vertex_base_2` is the default model for the concept 
`AlphaShapeVertex_2`. 

Parameters 
-------------- 

The template parameter `Traits` has to be a model of `AlphaShapeTraits_2`. 

The template parameter `Vb` has to be a model of `TriangulationVertexBase_2` (or `RegularTriangulationVertexBase_2`) 
if `Alpha_shape_vertex_base_2` is intended to be used with an alpha-shape class based on a 
`Delaunay_triangulation_2` (or a `Regular_triangulation_2`). 

The template parameter `ExactAlphaComparisonTag` is a tag that, when set to 
`CGAL::Tag_true`, triggers exact comparisons between alpha values. See the description 
provided in the documentation of `Alpha_shape_2` for more details. The default value is `CGAL::Tag_false`. 

\models ::AlphaShapeVertex_2 
*/
template< typename Traits, typename Vb, typename ExactAlphaComparisonTag >
class Alpha_shape_vertex_base_2 : public Vb {
public:

/// @}

}; /* end Alpha_shape_vertex_base_2 */
} /* end namespace CGAL */
