
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Triangulation_cell_base_3` is a model of the concept
`TriangulationCellBase_3`, the base cell of a 3D-triangulation.

This class can be used directly or can serve as a base to derive other classes
with some additional attributes (a color for example) tuned for a specific
application.


\tparam Traits is the geometric traits class and must be a model of `TriangulationTraits_3`.
It is actually not used by this class.

\tparam TDSCb is a combinatorial cell base class from which
`Triangulation_cell_base_3` derives.
It must be a model of the `TriangulationDSCellBase_3` concept.
It has the default value `Triangulation_ds_cell_base_3`.

\cgalModels `TriangulationCellBase_3`

\sa `CGAL::Triangulation_ds_cell_base_3`
\sa `CGAL::Triangulation_cell_base_with_info_3`
\sa `CGAL::Triangulation_vertex_base_3`

*/
template< typename Traits, typename TDSCb >
class Triangulation_cell_base_3 : public TDSCb {
public:

}; /* end Triangulation_cell_base_3 */
} /* end namespace CGAL */
