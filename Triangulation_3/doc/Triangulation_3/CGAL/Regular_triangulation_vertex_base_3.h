
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Regular_triangulation_vertex_base_3` is a model of the concept
`RegularTriangulationVertexBase_3`, the base vertex of a 3D regular triangulation.
This class stores a weighted point.

This class can be used directly or can serve as a base to derive other classes
with some additional attributes (a color for example) tuned for a specific
application.


\tparam Traits is the geometric traits class and must be a model of `RegularTriangulationTraits_3`.
Users of the geometric triangulations are strongly advised to use the same
geometric traits class as the one used in `Regular_triangulation_3`.
This way, the point type defined by the base vertex is
the same as the point type defined by the geometric traits class.

\tparam TDSVb is a combinatorial vertex base class from which
`Regular_triangulation_vertex_base_3` derives.
It must be a model of the `TriangulationDSVertexBase_3` concept.
It has the default value `Triangulation_ds_vertex_base_3<TDS>`.

\cgalModels{RegularTriangulationVertexBase_3}

\sa `CGAL::Triangulation_cell_base_3`
\sa `CGAL::Triangulation_ds_vertex_base_3`
\sa `CGAL::Triangulation_vertex_base_with_info_3`

*/
template< typename Traits, typename TDSVb >
class Regular_triangulation_vertex_base_3 : public TDSVb {
public:

/// \name Types
/// @{

/*!

*/
typedef Traits::Weighted_point_3 Point;

/// @}

}; /* end Regular_triangulation_vertex_base_3 */
} /* end namespace CGAL */
