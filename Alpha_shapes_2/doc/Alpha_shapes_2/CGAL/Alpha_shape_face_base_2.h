
namespace CGAL {

/*!
\ingroup PkgAlphaShape2

The class `Alpha_shape_face_base_2` is the default model for the concept `AlphaShapeFace_2`. 

Parameters 
-------------- 

The template parameter `Traits` has to be a model of `AlphaShapeTraits_2`. 

The template parameter `Fb` has to be a model of `TriangulationFaceBase_2` (or `RegularTriangulationFaceBase_2`) 
if `Alpha_shape_face_base_2` is intended to be used with an alpha-shape class based on a 
`Delaunay_triangulation_2` (or a `Regular_triangulation_2`). 

The template parameter `ExactAlphaComparisonTag` is a tag that, when set to 
`CGAL::Tag_true`, triggers exact comparisons between alpha values. See the description 
provided in the documentation of `Alpha_shape_2` for more details. The default value is `CGAL::Tag_false`. 

\models ::AlphaShapeFace_2 
*/
template< typename Traits, typename Fb, typename ExactAlphaComparisonTag >
class Alpha_shape_face_base_2 : public Fb {
public:

/// @}

}; /* end Alpha_shape_face_base_2 */
} /* end namespace CGAL */
