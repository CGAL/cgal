
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3Ref

The class `Alpha_shape_vertex_base_3` is the default model for the concept
`AlphaShapeVertex_3`.

\tparam Traits is the geometric traits class that is provided
to the `Alpha_shape_3` class.
\tparam Vb must be a vertex base class adapted to the type of triangulation that is being used.
        By default, it is instantiated with `Triangulation_vertex_base_3<Traits>`,
        which is appropriate for basic alpha shapes.
\tparam ExactAlphaComparisonTag is a tag that, when set to
\link Tag_true `Tag_true`\endlink, triggers exact comparisons between alpha values. See the description
provided in the documentation of `Alpha_shape_3` for more details. The default value is \link Tag_false `Tag_false`\endlink.
\tparam WeightedTag is used only if `ExactAlphaComparisonTag` is \link Tag_true `Tag_true`\endlink. It
must be \link Tag_true `Tag_true`\endlink if the underlying triangulation of the alpha shape to be used is a regular triangulation
and \link Tag_false `Tag_false`\endlink otherwise. The default is \link Tag_false `Tag_false`\endlink.

\cgalModels `AlphaShapeVertex_3`

\sa `Triangulation_vertex_base_3`
\sa `Regular_triangulation_vertex_base_3`
\sa `Periodic_3_triangulation_ds_vertex_base_3`
*/
template< typename Traits, typename Vb, typename ExactAlphaComparisonTag, typename WeightedTag >
class Alpha_shape_vertex_base_3 : public Vb {
public:


}; /* end Alpha_shape_vertex_base_3 */
} /* end namespace CGAL */
