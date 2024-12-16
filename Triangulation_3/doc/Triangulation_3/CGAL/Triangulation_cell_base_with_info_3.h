
namespace CGAL {

/*!
\ingroup PkgTriangulation3VertexCellClasses

The class `Triangulation_cell_base_with_info_3` is a model of the concept
`TriangulationCellBase_3`, the base cell of a 3D-triangulation.
It provides an easy way to add some user defined information in cells.
Note that input/output operators discard this additional information.


\tparam Info  is the information the user would like to add
to a cell. It has to be `DefaultConstructible` and `Assignable`.

\tparam TriangulationTraits_3  is the geometric traits class.
It is actually not used by this class.

\tparam TriangulationCellBase_3_ is a cell base class from which
`Triangulation_cell_base_with_info_3` derives.
It must be a model of the `TriangulationCellBase_3` concept.
It has the default value
`Triangulation_cell_base_3<TriangulationTraits_3>`.

\cgalModels{TriangulationCellBase_3,TriangulationCellBaseWithInfo_3}

\sa `CGAL::Triangulation_cell_base_3`
\sa `CGAL::Triangulation_vertex_base_with_info_3`

*/
template< typename Info, typename TriangulationTraits_3, typename TriangulationCellBase_3_ >
class Triangulation_cell_base_with_info_3 : public TriangulationCellBase_3_ {
public:

/// \name Types
/// @{

/*!

*/
typedef Info Info;

/// @}

/// \name Access Functions
/// @{

/*!
Returns a const reference to the object of type `Info` stored in the cell.
*/
const Info& info() const;

/*!
Returns a reference to the object of type `Info` stored in the cell.
*/
Info & info();

/// @}

}; /* end Triangulation_cell_base_with_info_3 */
} /* end namespace CGAL */
