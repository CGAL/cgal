
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

\cgalHeading{Parameters}

It is parameterized by base classes for vertices and cells which have to match 
the requirements for the concepts `TriangulationDSCellBase_3` and 
`TriangulationDSVertexBase_3` respectively. 

They have the default values `Triangulation_ds_vertex_base_3<TDS>` and 
`Triangulation_ds_cell_base_3<TDS>` respectively. 

The `Concurrency_tag` parameter allows to enable the use of a concurrent
container to store vertices and cells. It can be `Sequential_tag` (use of a 
`Compact_container` to store vertices and cells) or `Parallel_tag` 
(use of a `Concurrent_compact_container`). If it is 
`Parallel_tag`, the following functions can be called concurrently:
`create_vertex`, `create_cell`, `delete_vertex`, `delete_cell`.
`Sequential_tag` is the default value.

\cgalModels `TriangulationDataStructure_3`

The base class `Triangulation_utils_3` defines basic computations on
indices of vertices and neighbors of cells.

\attention All members listed here are additional to the interface
specified by the concept.

\sa `CGAL::Triangulation_ds_vertex_base_3` 
\sa `CGAL::Triangulation_ds_cell_base_3` 
\sa `CGAL::Triangulation_vertex_base_with_info_3` 
\sa `CGAL::Triangulation_cell_base_with_info_3` 
*/
template< typename TriangulationDSVertexBase_3, 
          typename TriangulationDSCellBase_3,
          typename Concurrency_tag >
class Triangulation_data_structure_3 : public CGAL::Triangulation_utils_3 {
public:

/// \name Types 
/// @{

/*!
Vertex container type. If Concurrency_tag is Parallel_tag, a
`Concurrent_compact_container` is used instead of a `Compact_container`.
*/ 
typedef Compact_container<Vertex, Default> Vertex_range; 

/*!
Cell container type. If Concurrency_tag is Parallel_tag, a
`Concurrent_compact_container` is used instead of a `Compact_container`.
*/ 
typedef Compact_container<Cell, Default> Cell_range; 
/// @} 

/// \name Operations 
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
