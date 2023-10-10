
namespace CGAL {

/*!
\ingroup PkgAlphaShapes3Ref

The class `Fixed_alpha_shape_cell_base_3` is the default model for the concept
`FixedAlphaShapeCell_3`.

\tparam Traits is the geometric traits class that is provided
to the `Alpha_shape_3` class.
\tparam Cb must be a cell base class adapted to the type of triangulation that is being used.
        By default, it is instantiated with `Delaunay_triangulation_cell_base_3<Traits>`,
        which is appropriate for basic alpha shapes.

\cgalModels{FixedAlphaShapeCell_3}

\sa `Alpha_shape_cell_base_3`
\sa `Delaunay_triangulation_cell_base_3`
\sa `Regular_triangulation_cell_base_3`
\sa `Periodic_3_triangulation_ds_cell_base_3`
*/
template< typename Traits, typename Cb >
class Fixed_alpha_shape_cell_base_3 : public Cb {
public:

}; /* end Fixed_alpha_shape_cell_base_3 */
} /* end namespace CGAL */
