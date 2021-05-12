
namespace CGAL {

/*!
\ingroup PkgTDS3Classes

The class `Triangulation_data_structure_3` stores a 3D-triangulation data
structure and provides the optional
geometric functionalities to be used as a parameter for a
3D-geometric triangulation (see Chapter \ref chapterTriangulation3 "3D Triangulations").

The vertices and cells are stored in two nested containers, which are
implemented using `Compact_container` (or `Concurrent_compact_container`,
see below). The class may offer some
flexibility for the choice of container in the future, in the form of
additional template parameters.

\tparam VertexBase must be a model of `TriangulationDSVertexBase_3`. The default is `Triangulation_ds_vertex_base_3<TDS>`.

\tparam CellBase  must be a model of `TriangulationDSCellBase_3`. The default is `Triangulation_ds_cell_base_3<TDS>`.

\tparam ConcurrencyTag enables the use of a concurrent
container to store vertices and cells. It can be `Sequential_tag` (use of a
`Compact_container` to store vertices and cells) or `Parallel_tag`
(use of a `Concurrent_compact_container`). If it is
`Parallel_tag`, the following functions can be called concurrently:
`create_vertex()`, `create_cell()`, `delete_vertex()`, and `delete_cell()`.
`Sequential_tag` is the default value.

\cgalModels `TriangulationDataStructure_3`

The base class `Triangulation_utils_3` defines basic computations on
indices of vertices and neighbors of cells.

\sa `CGAL::Triangulation_ds_vertex_base_3`
\sa `CGAL::Triangulation_ds_cell_base_3`
*/
template< typename VertexBase,
          typename CellBase,
          typename ConcurrencyTag >
class Triangulation_data_structure_3
  : public CGAL::Triangulation_utils_3
{
public:

/// \name Types
/// @{

typedef Triangulation_data_structure_2<VertexBase,FaceBase>  Tds;

/// The vertex type.
///
/// \sa Section \ref tds3cyclic
typedef  typename VertexBase::template Rebind_TDS<Tds>::Other  Vertex;

/// The face type.
///
/// \sa Section \ref tds3cyclic
typedef  typename CellBase::template Rebind_TDS<Tds>::Other  Cell;

/*!
Vertex container type. If `ConcurrencyTag` is `Parallel_tag`, a
`Concurrent_compact_container` is used instead of a `Compact_container`.
*/
typedef Compact_container<Vertex, Default> Vertex_range;

/*!
Cell container type. If `ConcurrencyTag` is `Parallel_tag`, a
`Concurrent_compact_container` is used instead of a `Compact_container`.
*/
typedef Compact_container<Cell, Default> Cell_range;
/// @}

/// \name Operations
///
/// In addition to the interface documented in the concept,
/// the class offers the following functions.
///
/// @{

/*!
Returns a reference to the container of cells.
*/
Cell_range & cells() const;

/*!
Returns a reference to the container of cells.
*/
Cell_range & cells();

/*!
Returns a reference to the container of vertices.
*/
Vertex_range & vertices() const;

/*!
Returns a reference to the container of vertices.
*/
Vertex_range & vertices();

/// @}

}; /* end Triangulation_data_structure_3 */
} /* end namespace CGAL */
