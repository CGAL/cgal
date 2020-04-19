
namespace CGAL {

/*!
\ingroup PkgAlphaShapes2Ref

The class `Alpha_shape_face_base_2` is the default model for the concept `AlphaShapeFace_2`.

\tparam Traits has to be a model of `AlphaShapeTraits_2`.

\tparam Fb has to be a model of `TriangulationFaceBase_2` (or `RegularTriangulationFaceBase_2`)
if `Alpha_shape_face_base_2` is intended to be used with an alpha-shape class based on a
`Delaunay_triangulation_2` (or a `Regular_triangulation_2`).

\tparam ExactAlphaComparisonTag is a tag that, when set to
\link Tag_true `Tag_true`\endlink, triggers exact comparisons between alpha values. See the description
provided in the documentation of `Alpha_shape_2` for more details. The default value is \link Tag_false `Tag_false`\endlink.

\cgalModels `AlphaShapeFace_2`

\sa `Triangulation_face_base_2`
\sa `Regular_triangulation_face_base_2`
\sa `Periodic_2_triangulation_face_base_2`
*/
template< typename Traits, typename Fb, typename ExactAlphaComparisonTag >
class Alpha_shape_face_base_2 : public Fb {
public:

}; /* end Alpha_shape_face_base_2 */
} /* end namespace CGAL */
