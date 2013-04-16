
namespace CGAL {

/*!
\ingroup PkgTDS3Classes

The class `Triangulation_data_structure_3` stores a 3D-triangulation data 
structure and provides the optional 
geometric functionalities to be used as a parameter for a 
3D-geometric triangulation (see Chapter \ref chapterTriangulation3 "3D Triangulations"). 

The vertices and cells are stored in two nested containers, which are 
implemented using `Compact_container`. The class may offer some 
flexibility for the choice of container in the future, in the form of 
additional template parameters. 

\cgalHeading{Parameters}

It is parameterized by base classes for vertices and cells which have to match 
the requirements for the concepts `TriangulationDSCellBase_3` and 
`TriangulationDSVertexBase_3` respectively. 

They have the default values `Triangulation_ds_vertex_base_3<>` and 
`Triangulation_ds_cell_base_3<>` respectively. 

The `Vertex_container_strategy` and `Cell_container_strategy` parameters are the 
strategies used by the `Compact_container` (or `Concurrent_compact_container` 
if the TDS is concurrency-safe) storing the vertices and cells. They are models of 
the concept `CompactContainerStrategy`. The default values are both
`Compact_container_strategy_base`.

The `Concurrency_tag` parameter allows to ask for a concurrency-safe TDS (with regard to
insertion and deletion of elements). Possible values are `CGAL::Sequential_tag` (the default) and
`CGAL::Parallel_tag`. The concurrency-safe version uses two `Concurrent_compact_container` to store 
vertices and cells (instead of two `Compact_container`).

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
          typename Vertex_container_strategy,
          typename Cell_container_strategy,
          typename Concurrency_tag >
class Triangulation_data_structure_3 : public CGAL::Triangulation_utils_3 {
public:

/// \name Types 
/// @{

/*! 
Vertex container type. If Concurrency_tag is Parallel_tag, a
`Concurrent_compact_container` is used instead of a `Compact_container`.
*/ 
typedef Compact_container<Vertex, Default, Vertex_container_strategy> Vertex_range; 

/*! 
Cell container type. If Concurrency_tag is Parallel_tag, a
`Concurrent_compact_container` is used instead of a `Compact_container`.
*/ 
typedef Compact_container<Cell, Default, Cell_container_strategy> Cell_range; 

/*! 
Strategy used by the `CGAL::Compact_container` (or `CGAL::Concurrent_compact_container` if the TDS
is concurrency-safe) storing the vertices.
*/ 
typedef Hidden_type Vertex_container_strategy;

/*! 
Strategy used by the `CGAL::Compact_container` (or `CGAL::Concurrent_compact_container` if the TDS
is concurrency-safe) storing the cells.
*/ 
typedef Hidden_type Cell_container_strategy;
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
