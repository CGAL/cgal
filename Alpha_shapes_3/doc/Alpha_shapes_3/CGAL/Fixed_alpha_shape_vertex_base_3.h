
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3Ref

The class `Fixed_alpha_shape_vertex_base_3` is the default model for the concept
`FixedAlphaShapeVertex_3`.

\tparam Traits is the geometric traits class that is provided
to the `Alpha_shape_3` class.
\tparam Vb must be a vertex base class adapted to the type of triangulation that is being used.
        By default, it is instantiated with `Triangulation_vertex_base_3<Traits>`,
        which is appropriate for basic alpha shapes.

\cgalModels{FixedAlphaShapeVertex_3}

\sa `Alpha_shape_vertex_base_3`
\sa `Triangulation_vertex_base_3`
\sa `Regular_triangulation_vertex_base_3`
\sa `Periodic_3_triangulation_ds_vertex_base_3`
*/
template< typename Traits, typename Vb >
class Fixed_alpha_shape_vertex_base_3 : public Vb {
public:

}; /* end Fixed_alpha_shape_vertex_base_3 */
} /* end namespace CGAL */
